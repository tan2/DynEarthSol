#include <iostream>
#include <limits>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "parameters.hpp"


static void declare_parameters(po::options_description &cfg,
                               Param &p)
{
    cfg.add_options()
        ("sim.modelname", po::value<std::string>(&p.sim.modelname),
         "Prefix for the output files")

        ("sim.max_steps", po::value<int>(&p.sim.max_steps),
         "Max. number of time steps")
        ("sim.max_time_in_yr", po::value<double>(&p.sim.max_time_in_yr),
         "Max. time (in years)")
        ("sim.output_step_interval", po::value<int>(&p.sim.output_step_interval),
         "Output step interval")
        ("sim.output_time_interval_in_yr", po::value<double>(&p.sim.output_time_interval_in_yr),
         "Output time interval (in years)")

        ("sim.is_restarting", po::value<bool>(&p.sim.is_restarting)->default_value(false),
         "Restarting from previous save?")
        ;

    cfg.add_options()
        ("mesh.meshing_option", po::value<int>(&p.mesh.meshing_option)->default_value(1),
         "How to create the new mesh?")

        ("mesh.xlength", po::value<double>(&p.mesh.xlength)->required(),
         "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(&p.mesh.ylength)->required(),
         "Length of y (in meters)")
        ("mesh.zlength", po::value<double>(&p.mesh.zlength)->required(),
         "Length of z (in meters)")

        ("mesh.resolution", po::value<double>(&p.mesh.resolution)->required(),
         "Spatial resolution (in meters)")
        // for 2D only
        ("mesh.min_angle", po::value<double>(&p.mesh.min_angle)->default_value(32.),
         "Min. angle of all triangles (in degrees)")
        // for 3D only
        ("mesh.min_tet_angle", po::value<double>(&p.mesh.min_tet_angle)->default_value(22.),
         "Min. dihedral angle of all tetrahedra (in degrees)")
        ("mesh.max_ratio", po::value<double>(&p.mesh.max_ratio)->default_value(2.),
         "Max. radius / length ratio of all tetrahedra")
        ;

    cfg.add_options()
        ("surface_temperature", po::value<double>(&p.surface_temperature)->default_value(273),
         "Surface temperature (in Kelvin)")

        ("mantle_temperature", po::value<double>(&p.mantle_temperature)->default_value(1600),
         "Mantle temperature (in Kelvin)")
        ;

}


static void read_parameters_from_file
(const char* filename,
 const po::options_description cfg,
 po::variables_map &vm)
{
    try {
        po::store(po::parse_config_file<char>(filename, cfg), vm);
        po::notify(vm);
    }
    catch (std::exception& e) {
        std::cerr << "Error reading config_file '" << filename << "'\n";
        std::cerr << e.what() << "\n";
        std::exit(1);
    }
}


static void validate_parameters(const po::variables_map &vm, Param &p)
{
    if ( ! (vm.count("sim.max_steps") || vm.count("sim.max_time_in_yr")) ) {
        std::cerr << "Must provide either sim.max_steps or sim.max_time_in_yr\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.max_steps") )
        p.sim.max_steps = std::numeric_limits<int>::max();;
    if ( ! vm.count("sim.max_time_in_yr") )
        p.sim.max_time_in_yr = std::numeric_limits<double>::max();;

    if ( ! (vm.count("sim.output_step_interval") || vm.count("sim.output_time_interval_in_yr")) ) {
        std::cerr << "Must provide either sim.output_step_interval or sim.output_time_interval_in_yr\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.output_step_interval") )
        p.sim.output_step_interval = std::numeric_limits<int>::max();;
    if ( ! vm.count("sim.output_time_interval_in_yr") )
        p.sim.output_time_interval_in_yr = std::numeric_limits<double>::max();;
}


void get_input_parameters(const char* filename, Param& p)
{
    po::options_description cfg("Config file options");
    po::variables_map vm;

    declare_parameters(cfg, p);
    read_parameters_from_file(filename, cfg, vm);
    validate_parameters(vm, p);
}
