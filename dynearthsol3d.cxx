
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <limits>
#include <iostream>
#include <fstream>
//#include <string>

#include "parameters.hpp"

void get_input_parameters(const char* filename, Param& param)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "Error: cannot open config_file '" << filename << "'\n";
        std::exit(1);
    }

    //
    // declare input parameters
    //
    po::options_description cfg("Config file options");
    cfg.add_options()
        ("sim.max_steps", po::value<int>(&param.sim.max_steps), "Max. number of time steps")
        ("sim.max_time", po::value<double>(&param.sim.max_time), "Max. time (in seconds)")
        ("sim.output_step_interval", po::value<int>(&param.sim.output_step_interval),
         "Output step interval")
        ("sim.output_time_interval", po::value<double>(&param.sim.output_time_interval),
         "Output time interval")
        ("sim.is_restarting", po::value<bool>(&param.sim.is_restarting)->default_value(false),
         "Restarting from previous save?")
        ;

    cfg.add_options()
        ("mesh.xlength", po::value<double>(&param.mesh.xlength)->required(), "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(&param.mesh.ylength)->required(), "Length of y (in meters)")
        ("mesh.zlength", po::value<double>(&param.mesh.zlength)->required(), "Length of z (in meters)")
        ;

    //
    // parse config file
    //
    po::variables_map vm;
    try {
        po::store(po::parse_config_file(ifs, cfg), vm);
        po::notify(vm);
    }
    catch (std::exception& e) {
        std::cerr << "Error reading config_file '" << filename << "'\n";
        std::cerr << e.what() << "\n";
        std::exit(1);
    }

    //
    // validate parameters
    //
    if ( ! (vm.count("sim.max_steps") || vm.count("sim.max_time")) ) {
        std::cerr << "Must provide either sim.max_steps or sim.max_time\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.max_steps") )
        param.sim.max_steps = std::numeric_limits<int>::max();;
    if ( ! vm.count("sim.max_time") )
        param.sim.max_time = std::numeric_limits<double>::max();;

    if ( ! (vm.count("sim.output_step_interval") || vm.count("sim.output_time_interval")) ) {
        std::cerr << "Must provide either sim.output_step_interval or sim.output_time_interval\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.output_step_interval") )
        param.sim.output_step_interval = std::numeric_limits<int>::max();;
    if ( ! vm.count("sim.output_time_interval") )
        param.sim.output_time_interval = std::numeric_limits<double>::max();;

    return;
}

void init() {};
void restart() {};
void update_temperature() {};
void update_strain_rate() {};
void update_stress() {};
void update_force() {};
void update_mesh() {};
void rotate_stress() {};
void output() {};


int main(int argc, const char* argv[])
{
    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        return -1;
    }

    Param param;
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    int steps = 0;
    double time = 0;
    if (! param.sim.is_restarting) {
        init();
        output();
    }
    else {
        restart();
    }

    int frame = 1;
    do {
        double dt = 1e7;
        steps++;
        time += dt;

        update_temperature();
        update_strain_rate();
        update_stress();
        update_force();
        update_mesh();
        rotate_stress();

        std::cout << "Step: " << steps << ", time:" << time << "\n";

        if ( (steps >= frame*param.sim.output_step_interval) ||
             (time >= frame*param.sim.output_time_interval) ) {
            output();
            std::cout << frame <<"-th output\n";
            frame++;
        }

    } while (steps < param.sim.max_steps && time <= param.sim.max_time);

    return 0;
}
