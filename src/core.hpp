/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

//vtk includes
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>

//std includes

#include <algorithm>
#include <chrono>
#include <chrono>
#include <cstdlib>
#include <cstdlib>
#include <errno.h>
#include <fstream>
#include <map>
#include <memory> //unique ptr
#include <set>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h> //for getpid
#include <utility> // std::pair
#include <vector>

//boost includes
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/bind/bind.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
namespace pt = boost::property_tree;
namespace po = boost::program_options;

// tbb
#include <tbb/concurrent_vector.h>

//osgeo
#include <ogr_spatialref.h>

//gls
#include <gsl/gsl_errno.h>

//includes from CHM
#include "exception.hpp"
#include "filter_base.hpp"
#include "global.hpp"
#include "interpolation.hpp"
#include "logger.hpp"
#include "math/coordinates.hpp"
#include "metdata.hpp"
#include "module_base.hpp"
#include "readjson.hpp"
#include "station.hpp"
#include "str_format.h"
#include "timer.hpp"
#include "timeseries/netcdf.hpp"
#include "triangulation.hpp"
#include "version.h"

#ifdef USE_MPI
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#endif

struct vertex{
    std::string name;
};

struct edge{
    std::string variable;
};

//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
//      boost::property<boost::vertex_color_t, boost::default_color_type>
//    > Graph;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,boost::property<boost::vertex_index_t,int,vertex>,edge> Graph;


//from http://stackoverflow.com/questions/11369115/how-to-print-a-graph-in-graphviz-with-multiple-properties-displayed
template <class VariableMap>
class edge_writer {
public:
    edge_writer(VariableMap v) : vm(v) {}
    template <class Edge>
    void operator()(ostream &out, const Edge& e) const {
        out << "[label=\"" << vm[e] << "\", edgetype=" << vm[e] << "]";
    }
private:
    VariableMap vm;
};

template <class VariableMap>
inline edge_writer<VariableMap>
make_edge_writer(VariableMap v) {
    return edge_writer<VariableMap>(v);
}

//typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
//typedef std::pair<int, int> Edge;


/**
 * The main model core, handles initialization of the model
**/

class core
{
    friend class CoreTest;
public:


  /**
   * Reads the main JSON configuration file. It assumes the base of the JSON is an object. That is, the file
   * starts with { ... }.
   * Within this file are a collection of meshes that are expected to have the same number of x,y
   * points. This is do
   *
   * ne so that, for example, elevation, forest cover, sky-view factor, etc
   * may be added individually. Generation of the meshes should be done via the utilities for this.
   * An example of mesh.config is:
   * \code
   *  {
   *    "meshes":
   *    {
   *            "DEM":
   *            {
   *                    "file": "mesh.asc"
   *            }
   *            ,
   *            "Veg":
   *            {
   *                    "file": "veg.asc"
   *            },
   *            "svf":
   *            {
   *                    "file": "svf.asc"
   *            }
   *    }
   *   }
   *   \endcode
   * @param file The file to open
  **/
    void init(int argc, char **argv);
    void config_options( pt::ptree &value);
    void config_modules( pt::ptree& value,const pt::ptree& config,std::vector<std::string> remove,std::vector<std::string> add);

    /**
     *
     * @param value True if loading a partitioned mesh
     * @return
     */
    bool config_meshes( pt::ptree& value);
    void config_forcing(pt::ptree& value);
    void config_module_overrides( pt::ptree& value);
    void config_parameters(pt::ptree &value);
    void config_output(pt::ptree& value);
    void config_global( pt::ptree& value);
    void config_checkpoint( pt::ptree& value);

    /**
     * Determines what the start end times should be, and ensures consistency from a check pointed file
     */
    void determine_startend_ts_forcing();

    /**
     * Determines the order modules need to be scheduleled in to maximize parallelism
     */
    void _schedule_modules();
    void _find_and_insert_subjson(pt::ptree& value);

    /**
     * Populates a list of stations needed within each face
     */
    void populate_face_station_lists();

    /**
     * Populates a list of stations needed on each MPI process
     */
    void populate_distributed_station_lists();

    /**
     * Checks if the mesh is geographic
     * @param path
     */
    bool check_is_geographic(const std::string& path);

    // .first = config file to use
    // .second = extra options, if any.
    typedef boost::tuple<
            std::string, //config file path to load. defaults to CHM.config
            std::vector<std::pair<std::string,std::string>>, //insert or overide config value
            std::vector<std::string>, //remove config value
            std::vector<std::string>, //remove module
            std::vector<std::string>,  // add module
            bool //legacy-log
    > cmdl_opt;

    cmdl_opt config_cmdl_options(int argc, char **argv);

    /**
     * Initializes the logger
     */
    core();
    ~core();

    void run();

    /**
     * Shutdown. In MPI mode allows us to trigger an MPI_Abort on exception
     * @param abort
     */
    void end(const bool abort = false );

    std::vector< std::pair<module,size_t> >& get_active_module_list();

    pt::ptree _cfg;
    boost::filesystem::path output_folder_path; //path to output folder
    boost::filesystem::path log_file_path; // fully qualified path to the log file

    // On a clean exit, a sentinel file will be written
    bool clean_exit;
protected:

    std::string version = "CHM " GIT_BRANCH "/" GIT_COMMIT_HASH;

    //current level of the logger. Defaults to debug, but can be changed via configuration settings
    log_level _log_level;

    // if the users passes in a config file path that isn't the current directory
    // e.g. CHM -f /some/other/path/CHM.json
    // then we need to affix every file IO (excep the log ?) with this path.
    boost::filesystem::path cwd_dir;


    bool _output_station_ptv; //should we output the station ptv file? if we have no output section, then don't do this.

    //this is called via system call when the model is done to notify the user
    std::string _notification_script;

    //main mesh object
    boost::shared_ptr< triangulation > _mesh;

    // these are saved here so-as to be used elsewhere
    std::string _mesh_path;


    //if radius selection for stations is chosen this holds that
    double radius;
    double N; // meters, radius for station search

    //holds all the modules that are to be run on each mesh element
    //pair as we also need to store the make order
    std::vector< std::pair<module,size_t> > _modules;
    std::vector< std::vector < module> > _chunked_modules;
    std::vector< std::pair<std::string,std::string> > _overrides;
    boost::shared_ptr<global> _global;

    bool _use_netcdf; // flag if we are using netcdf. If we are, it enables incremental reads of the netcdf file for speed.
    std::shared_ptr<metdata> _metdata; //met data loader, shared for use with boost::bind

    //calculates the order modules are to be run in
    void _determine_module_dep();

    interp_alg _interpolation_method;

    //holds a unique list of all variables provided by all the met files;
    std::set<std::string> _provided_var_met_files;
    //unique list of all variables provided by all the modules
    std::set<std::string> _provided_var_module;
    std::set<std::string> _provided_var_vector;

    //unique set of all the paramters provided by the meshes
    std::set<std::string> _provided_parameters;
    std::set<std::string> _provided_initial_conditions;


    boost::posix_time::ptime* _start_ts;
    boost::posix_time::ptime* _end_ts;

    struct point_mode_info
    {
        bool enable;

        // The default mode of point mode is to use whatever stations we'd use for the face containing
        //  this output. If we ask sepficially for a single station, then only that station will be used.
        bool use_specific_station;
        std::string forcing; // empty unless the above is set true

    } point_mode;


    class output_info
    {
    public:
        output_info():
                      fname{""},
                      latitude{0}, longitude{0},
                      name{""},
                      x{0}, y{0},
                      only_last_n{SIZE_MAX}
        {
            face = nullptr;
        }

        enum output_type
        {
            time_series,
            mesh
        };
        enum mesh_outputs
        {
            vtp,
            vtu,
            ascii
        };

        // Should we output?
        bool should_output(const size_t& max_ts,
                           const size_t& current_ts,
                           const boost::posix_time::ptime& _current_date
                           )
        {
            bool should_output = false;

            if(only_last_n)
            {
                auto ts_left = max_ts - current_ts;
                if( ts_left <= *only_last_n) // if we are within the last n timesteps, output
                    should_output = true;
            }

            if(frequency)
            {
                if(current_ts % *frequency == 0)
                    should_output = true;
            }

            if(specific_time)
            {
                if( (_current_date.time_of_day().hours() == specific_time->time_of_day().hours()) &&
                    (_current_date.time_of_day().minutes() == specific_time->time_of_day().minutes()) )
                    should_output = true;
            }

            if(specific_datetime)
            {
                if(_current_date == *specific_datetime)
                    should_output = true;
            }

            return should_output;

        }

        // print to stdout DEBUG all the valid outputs selected
        void list_outputs()
        {
                SPDLOG_DEBUG("Output frequency options for {}", name);

                if(only_last_n)
                    SPDLOG_DEBUG("\tonly_last_n = {}", *only_last_n);
                if(frequency)
                    SPDLOG_DEBUG("\tfrequency = {}", *frequency);
                if(specific_time)
                    SPDLOG_DEBUG("\tspecific_time = {}", std::to_string(specific_time->time_of_day().hours()) + ":" + std::to_string(specific_time->time_of_day().minutes()));
                if(specific_datetime)
                    SPDLOG_DEBUG("\tspecific_datetime = {}", boost::posix_time::to_simple_string(*specific_datetime));
        }

        output_type type; // the type of output
        std::string name;
        std::vector<mesh_outputs> mesh_output_formats;
        std::string fname;

        // these are input by the user, assumed to be WGS84
        double latitude;
        double longitude;

        // if we are outputting on a projected mesh then we need to store the projected coords here
        double x;
        double y;

        std::set<std::string> variables;
        mesh_elem face;
        timeseries ts;

        // Output options

        // every n timesteps
        boost::optional<size_t> frequency;

        // at a specific date-time
        boost::optional<boost::posix_time::ptime> specific_datetime;

        // at a specific time
        boost::optional<boost::posix_time::ptime> specific_time;


        //Only output the last n timesteps. -1 = all
        boost::optional<size_t> only_last_n;

    };

    std::vector<output_info> _outputs;


    // Detects various information about the HPC scheduler we might be run der
    class hpc_scheduler_info
    {
    public:

      boost::posix_time::time_duration max_wallclock; // maximum wallclock in seconds
      boost::posix_time::ptime wallclock_start; // time we started the simulation at
      bool has_wallclock_limit;

        hpc_scheduler_info()
        {
            max_wallclock = boost::posix_time::seconds(0);
            has_wallclock_limit = false;
        }

        /**
         * If we have a wallclock limit, how much time left?
         * Only produces a useful delta if has_wallclock_limit = true;
         * @return
         */
        boost::posix_time::time_duration wallclock_remaining()
        {
            return max_wallclock - (boost::posix_time::second_clock::local_time()-wallclock_start);
        }

        void detect()
        {
            // Check if we are running under slurm
            const char* SLURM_JOB_ID = std::getenv("SLURM_JOB_ID");
            if (SLURM_JOB_ID)
            {
                const char* SLURM_TASK_PID = std::getenv("SLURM_TASK_PID"); // The process ID of the task being started.
                const char* SLURM_PROCID =
                    std::getenv("SLURM_PROCID"); // The MPI rank (or relative process ID) of the current process

                SPDLOG_DEBUG("Detected running under SLURM as jobid {}", SLURM_JOB_ID);
                SPDLOG_DEBUG("SLURM_TASK_PID = {}", SLURM_TASK_PID);
                SPDLOG_DEBUG("SLURM_PROCID = {} ", SLURM_PROCID);
            }


            // check if we are running under PBS
            const char* PBS_JOB_ID = std::getenv("PBS_JOBID");
            if(PBS_JOB_ID)
            {
                SPDLOG_DEBUG("Detected running under PBS as jobid {}", PBS_JOB_ID);
            }

            const char* CHM_WALLCLOCK = std::getenv("CHM_WALLCLOCK_LIMIT");
            if(CHM_WALLCLOCK)
            {
                try {
                    max_wallclock = boost::posix_time::duration_from_string(CHM_WALLCLOCK);
                    has_wallclock_limit = true;
                    wallclock_start = boost::posix_time::second_clock::local_time();
                    SPDLOG_DEBUG("Detected a max wallclock of {}", boost::posix_time::to_simple_string(max_wallclock));
                } catch (...) {
                    CHM_THROW_EXCEPTION(chm_error, "The value given for environment variable CHM_WALLCLOCK is invalid");
                }
            }

        }
    } _hpc_scheduler_info;


    // Checkpointing options
    class chkptOp
    {
      public:
        chkptOp():
                    do_checkpoint{false},
                    load_from_checkpoint{false},
                    on_last{false},
                    checkpoint_request_terminate{false}
        {
            abort_when_wallclock_left = boost::posix_time::minutes(2);
        }

        boost::filesystem::path ckpt_path; // root path to chckpoint folder
        netcdf in_savestate; // if we are loading from checkpoint
        bool do_checkpoint; // should we check point?
        bool load_from_checkpoint; // are we loading from a checkpoint?
        // amount of time to give ourselves to bail and checkpoint if we have a wall clock limit
        boost::posix_time::time_duration abort_when_wallclock_left;
        boost::optional<bool> on_outta_time; // bail when we are out of time
        boost::optional<bool> on_last; //only checkpoint on the last timestep
        boost::optional<size_t> frequency; // frequency of checkpoints


        // used to stop the simulation when we checkpoint when we are outta time
        bool checkpoint_request_terminate;

        /**
         * Should checkpointing occur
         * @param current_ts
         * @param is_last_ts
         * @return
         */
        bool should_checkpoint(size_t current_ts, bool is_last_ts, hpc_scheduler_info& scheduler_info, boost::mpi::communicator& comm_world)
        {
            if(!do_checkpoint)
                return false;

            if(on_last && *on_last && is_last_ts)
                return true;

            // don't checkpoint on the first ts if we are doing frequency checkpoints
            if( frequency && current_ts !=0 && (current_ts % *frequency ==0) )
                return true;

            // check if we are running out of time
            if(on_outta_time && *on_outta_time &&
                scheduler_info.has_wallclock_limit 
                )
            {
                
                int outoftime =  scheduler_info.wallclock_remaining() <= abort_when_wallclock_left;
                int global_outoftime = -1;

                // find is anyone thinks we should bail
                boost::mpi::all_reduce(comm_world, outoftime, global_outoftime, boost::mpi::maximum<int>());

                if(global_outoftime)
                {
                    SPDLOG_DEBUG("Detected wallclock of {} remaining. Triggering checkpoint.",
                                boost::posix_time::to_simple_string(scheduler_info.wallclock_remaining()));
                    checkpoint_request_terminate = true;
                    return true;
                }

            }

            return false;
        }


    } _checkpoint_opts;



    //command line argument options we need to keep track of

    struct
    {
        bool tmp;  // empty until we use this more

    } cli_options;


#ifdef USE_MPI
    boost::mpi::environment _mpi_env;
    boost::mpi::communicator _comm_world;
#endif

};

// Macro to quickyl check return of system calls
// TODO: Decide what to do if system call error has occurred
#define CHK_SYSTEM_ERR(ierr)					\
  if (ierr < 0) {						\
    SPDLOG_ERROR(strerror(errno));				\
  };
