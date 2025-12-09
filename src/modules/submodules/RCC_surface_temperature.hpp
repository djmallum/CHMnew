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
#include <algorithm>
#include <concepts>

/**
 * @brief Radiative-Conductive-Convective (RCC) ground surface temperature model.
 *
 * This submodule calculates ground surface temperature using the RCC approach,
 * which employs an empirical formula derived from linear regression analysis.
 * The model combines air temperature, net radiation, and thaw front depth to
 * estimate surface temperature, particularly useful for modeling frost table
 * dynamics during thaw seasons.
 *
 * The implementation is based on Equation (7) from:
 * T. J. Williams, J. W. Pomeroy, J. R. Janowicz, S. K. Carey, K. Rasouli,
 * and W. L. Quinton, "A radiative–conductive–convective approach to calculate
 * thaw season ground surface temperatures for modelling frost table dynamics,"
 * Hydrological Processes, vol. 29, no. 18, pp. 3954–3965, Aug. 2015.
 *
 * Key features:
 * - Uses arctan function for smooth transition between temperature regimes
 * - Tracks maximum thaw front depth over time
 * - Combines radiative, conductive, and convective heat transfer components
 *
 * Usage example:
 * @code
 * struct MyData {
 *     // Implement all RCC_data concept requirements...
 * };
 *
 * RCC_surface_temperature<MyData> rcc_model;
 * MyData data;
 * rcc_model.execute(data);  // Updates data.surface_temperature()
 * @endcode
 */

template<typename T>
concept RCC_data = requires(T& t)
{
    // Inputs
	{ t.thaw_front_depth() } -> std::floating_point;

	{ t.air_temperature() } -> std::floating_point; 

	{ t.net_radiation() } -> std::floating_point; 
    
    { t.thaw_front_depth_last() } -> std::floating_point;
    
    // Outputs    
	{ t.surface_temperature(std::declval<const double>()) } -> std::same_as<void>;

    { t.thaw_front_depth_last(std::declval<const double>()) } -> std::same_as<void>;
};

template<RCC_data data>
class RCC_surface_temperature : public base_step<RCC_surface_temperature<data>,data>
{
public:
	explicit RCC_surface_temperature() {};
	~RCC_surface_temperature() {};

	void execute(data& d) override final;

private:
    static constexpr auto _A = 0.77;
    static constexpr auto _B = 0.02;
    static constexpr auto _C = 7.0;
    static constexpr auto _D = 0.03; // d is taken by input argument
    static constexpr auto _ARCTAN_NORMALIZER = 2.0 / 3.14156265; 

};


template<RCC_data data>
void RCC_surface_temperature<data>::execute(data& d)
{
    // arctan is often used because it varies smoothly from -pi/2 to pi/2
    // this value makes it from -1 to 1.
    // TODO Possible bug: thaw_front_depth_last is never reset
    // This bug exists in CRHM too
    d.thaw_front_depth_last(std::max(d.thaw_front_depth_last(),
            d.thaw_front_depth()));
    

	auto T = 
		(_A * d.air_temperature() + _B * d.net_radiation()) * 
		std::atan(_C * (d.thaw_front_depth_last() + _D)) * _ARCTAN_NORMALIZER;
	d.surface_temperature(T);
};

