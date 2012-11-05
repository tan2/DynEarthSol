#include <cstdio>
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

        // for meshing_option = 2 only
        /* read the value as string, then parse as two numbers later */
        ("mesh.refined_zonex", po::value<std::string>(),
         "Refining portion of xlength (two numbers between 0~1)")
        ("mesh.refined_zoney", po::value<std::string>(),
         "Refining portion of ylength (two numbers between 0~1)")
        ("mesh.refined_zonez", po::value<std::string>(),
         "Refining portion of zlength (two numbers between 0~1)")
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

    if (p.mesh.meshing_option == 2) {
        if ( ! vm.count("mesh.refined_zonex") ||
#ifdef THREED
             ! vm.count("mesh.refined_zoney") ||
#endif
             ! vm.count("mesh.refined_zonez") ) {
        std::cerr << "Must provide mesh.refined_zonex, "
#ifdef THREED
                  << "mesh.refined_zoney, "
#endif
                  << "mesh.refined_zonez.\n";
        std::exit(1);
        }

        /* get 2 numbers from the string */
        double d0, d1;
        int n;
        std::string str;
        str = vm["mesh.refined_zonex"].as<std::string>();
        n = std::sscanf(str.c_str(), "[%lf, %lf]", &d0, &d1);
        if (n != 2 || d0 < 0 || d0 > 1 || d1 < 0 || d1 > 1) {
            std::cerr << "Error: incorrect value for mesh.refine_zonex,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0,d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonex.first = d0;
        p.mesh.refined_zonex.second = d1;
#ifdef THREED
        str = vm["mesh.refined_zoney"].as<std::string>();
        n = std::sscanf(str.c_str(), "[%lf, %lf]", &d0, &d1);
        if (n != 2 || d0 < 0 || d0 > 1 || d1 < 0 || d1 > 1) {
            std::cerr << "Error: incorrect value for mesh.refine_zoney,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0,d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zoney.first = d0;
        p.mesh.refined_zoney.second = d1;
#endif
        str = vm["mesh.refined_zonez"].as<std::string>();
        n = std::sscanf(str.c_str(), "[%lf, %lf]", &d0, &d1);
        if (n != 2 || d0 < 0 || d0 > 1 || d1 < 0 || d1 > 1) {
            std::cerr << "Error: incorrect value for mesh.refine_zonez,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0,d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonez.first = d0;
        p.mesh.refined_zonez.second = d1;
    }
}


void get_input_parameters(const char* filename, Param& p)
{
    po::options_description cfg("Config file options");
    po::variables_map vm;

    declare_parameters(cfg, p);
    read_parameters_from_file(filename, cfg, vm);
    validate_parameters(vm, p);
}
