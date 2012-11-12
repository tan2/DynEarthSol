#include <cstdio>
#include <iostream>
#include <limits>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "parameters.hpp"
#include "matprops.hpp"


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
         "Length of y (in meters), for 3D only")
        ("mesh.zlength", po::value<double>(&p.mesh.zlength)->required(),
         "Length of z (in meters)")

        ("mesh.resolution", po::value<double>(&p.mesh.resolution)->required(),
         "Spatial resolution (in meters)")
        // for 2D only
        ("mesh.min_angle", po::value<double>(&p.mesh.min_angle)->default_value(32.),
         "Min. angle of all triangles (in degrees), for 2D only")
        // for 3D only
        ("mesh.min_tet_angle", po::value<double>(&p.mesh.min_tet_angle)->default_value(22.),
         "Min. dihedral angle of all tetrahedra (in degrees), for 3D only")
        ("mesh.max_ratio", po::value<double>(&p.mesh.max_ratio)->default_value(2.),
         "Max. radius / length ratio of all tetrahedra, for 3D only")

        // for meshing_option = 2 only
        /* read the value as string, then parse as two numbers later */
        ("mesh.refined_zonex", po::value<std::string>(),
         "Refining portion of xlength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only")
        ("mesh.refined_zoney", po::value<std::string>(),
         "Refining portion of ylength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only, for 3D only")
        ("mesh.refined_zonez", po::value<std::string>(),
         "Refining portion of zlength ([d0,d1]; 0<=d0<=d1<=1), for meshing_option=2 only")
        ;

    cfg.add_options()
        ("control.gravity", po::value<double>(&p.control.gravity)->default_value(10),
         "Magnitude of the gravity (in m/s^2)")

        ("control.inertial_scaling", po::value<double>(&p.control.inertial_scaling)->default_value(1e5),
         "Scaling factor for inertial (a large number)")
        ;

    cfg.add_options()
        ("bc.surface_temperature", po::value<double>(&p.bc.surface_temperature)->default_value(273),
         "Surface temperature (in Kelvin)")
        ("bc.mantle_temperature", po::value<double>(&p.bc.mantle_temperature)->default_value(1600),
         "Mantle temperature (in Kelvin)")
        ("bc.max_vbc_val", po::value<double>(&p.bc.max_vbc_val)->default_value(1e-9),
         "Magnitude of boundary velocity (in m/s)")
        ;

    cfg.add_options()
        ("mat.rheology_type", po::value<std::string>()->required(),
         "Type of rheology, either 'elastic', 'viscous', 'maxwell', 'elasto-plastic', or 'elasto-viscous-plastic'")
        ("mat.num_material", po::value<int>(&p.mat.nmat)->default_value(1),
         "Number of material types")
        ("mat.max_viscosity", po::value<double>(&p.mat.visc_max)->default_value(1e24),
         "Max. value of viscosity")
        ("mat.min_viscosity", po::value<double>(&p.mat.visc_min)->default_value(1e18),
         "Min. value of viscosity")
        ("mat.max_thermal_diffusivity", po::value<double>(&p.mat.therm_diff_max)->default_value(5e-6),
         "Max. value of thermal diffusivity")
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
    //
    // stopping condition and output interval are based on either model time or step
    //
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


    //
    // these parameters are required in mesh.meshing_option == 2
    //
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
        if (n != 2 || d0 < 0 || d1 > 1 || d0 > d1) {
            std::cerr << "Error: incorrect value for mesh.refine_zonex,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonex.first = d0;
        p.mesh.refined_zonex.second = d1;
#ifdef THREED
        str = vm["mesh.refined_zoney"].as<std::string>();
        n = std::sscanf(str.c_str(), "[%lf, %lf]", &d0, &d1);
        if (n != 2 || d0 < 0 || d1 > 1 || d0 > d1) {
            std::cerr << "Error: incorrect value for mesh.refine_zoney,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zoney.first = d0;
        p.mesh.refined_zoney.second = d1;
#endif
        str = vm["mesh.refined_zonez"].as<std::string>();
        n = std::sscanf(str.c_str(), "[%lf, %lf]", &d0, &d1);
        if (n != 2 || d0 < 0 || d1 > 1 || d0 > d1) {
            std::cerr << "Error: incorrect value for mesh.refine_zonez,\n"
                      << "       must in this format '[d0, d1]', 0 <= d0 <= d1 <= 1.\n";
            std::exit(1);
        }
        p.mesh.refined_zonez.first = d0;
        p.mesh.refined_zonez.second = d1;
    }

    //
    // material properties
    //
    {
        std::string str = vm["mat.rheology_type"].as<std::string>();
        if (str == std::string("elastic"))
            p.mat.rheol_type = MatProps::rh_elastic;
        else if (str == std::string("viscous"))
            p.mat.rheol_type = MatProps::rh_viscous;
        else if (str == std::string("maxwell"))
            p.mat.rheol_type = MatProps::rh_maxwell;
        else if (str == std::string("elasto-plastic"))
            p.mat.rheol_type = MatProps::rh_ep;
        else if (str == std::string("elasto-viscous-plastic"))
            p.mat.rheol_type = MatProps::rh_evp;
        else {
            std::cerr << "Error: unknown rheology: '" << str << "'\n";
            std::exit(1);
        }
    }

}


void get_input_parameters(const char* filename, Param& p)
{
    po::options_description cfg("Config file options");
    po::variables_map vm;

    declare_parameters(cfg, p);
    // print help message
    if (std::strncmp(filename, "-h", 3) == 0 ||
        std::strncmp(filename, "--help", 7) == 0) {
        std::cout << cfg;
        std::exit(0);
    }
    read_parameters_from_file(filename, cfg, vm);
    validate_parameters(vm, p);
}
