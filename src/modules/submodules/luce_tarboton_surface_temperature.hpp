//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//
#pragma once
#include "base_step.hpp"
#include <cmath>
#include <concepts>
/**
 * @brief Luce-Tarboton snow surface temperature model.
 *
 * This submodule calculates surface temperature below the snow surface using 
 * the Luce-Tarboton approach, which solves the heat conduction equation through 
 * snowpack. The model accounts for snow density effects on thermal conductivity 
 * and switches between low-density and high-density thermal conductivity formulas.
 *
 * Based on Equation (11) from:
 * C. H. Luce and D. G. Tarboton, "Evaluation of alternative formulae for 
 * calculation of surface temperature in snowmelt models using frequency 
 * analysis of temperature observations," Hydrol. Earth Syst. Sci., 2010.
 *
 * Note: This implementation follows CRHM's adaptation, using snow_depth/2 
 * instead of the damping depth parameter 'd' from the original paper.
 *
 * Key features:
 * - Thermal conductivity models based on snow density threshold (156 kg/m³)
 * - Fallback to air temperature for extreme cold conditions (< -70°C mean temp)
 * - Unit conversions between kg/m³ and g/cm³ for thermal conductivity formulas
 *
 * Usage example:
 * @code
 * struct MyData {
 *     // Implement all luce_tarboton_data concept requirements...
 * };
 *
 * luce_tarboton_surface_temperature<MyData> lt_model;
 * MyData data;
 * lt_model.execute(data);  // Updates surface_temperature() and snow_thermal_conductivity()
 * @endcode
 */
  
template<class T>
concept luce_tarboton_data = requires(T& t)
{
    // Inputs	
    { t.air_temperature() } -> std::floating_point;

	{ t.snow_depth() } -> std::floating_point;

	{ t.snow_density() } -> std::floating_point;

	{ t.ground_heat_flux() } -> std::floating_point;

    { t.daily_mean_temperature() } -> std::floating_point;

    // Outputs
	{ t.surface_temperature(std::declval<const double>())} -> std::same_as<void>;

	{ t.snow_thermal_conductivity(std::declval<const double>()) } -> std::same_as<void>;
};

template<luce_tarboton_data data>
class luce_tarboton_surface_temperature : public base_step<luce_tarboton_surface_temperature<data>,data>
{
public:
	explicit luce_tarboton_surface_temperature() {};
	~luce_tarboton_surface_temperature() {};

	void execute_impl(data& d);

private:

    static constexpr auto kg_per_m3_to_g_per_cm3 = 1.0 / 1000.0;
    static constexpr auto density_threshold = 156.0;
    static constexpr auto temperature_threshold = -70.0;

    double low_density_thermal_conductivity(data& d);  
    double high_density_thermal_conductivity(data& d);
};

template<luce_tarboton_data data>
void luce_tarboton_surface_temperature<data>::execute_impl(data& d)
{
    
	double tc,T;
	if (d.snow_density() < density_threshold)
	{
        tc = low_density_thermal_conductivity(d);
	}
	else
	{
		tc = high_density_thermal_conductivity(d); 
	}

	d.snow_thermal_conductivity(tc);

	if (d.daily_mean_temperature() < temperature_threshold)
		T = d.air_temperature();
	else
	{
		T = d.daily_mean_temperature() + 
			d.ground_heat_flux() * 0.5 * d.snow_depth() / tc;
	}

	d.surface_temperature(T);
};

template<luce_tarboton_data data>
double luce_tarboton_surface_temperature<data>::low_density_thermal_conductivity(data& d)
{
    constexpr auto a = 0.023;
    constexpr auto b = 0.234;

    return a + b * d.snow_density() * kg_per_m3_to_g_per_cm3;
};

template<luce_tarboton_data data>
double luce_tarboton_surface_temperature<data>::high_density_thermal_conductivity(data& d)
{
    constexpr auto a = 0.138;
    constexpr auto b = 1.01;
    constexpr auto c = 3.233;

    double snow_density = d.snow_density() * kg_per_m3_to_g_per_cm3;

    return a - b * snow_density + c * std::pow(snow_density,2);
};

