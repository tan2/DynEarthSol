
#include <boost/program_options.hpp>
namespace po = boost::program_options;

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
    if ( !(vm.count("sim.max_steps") || vm.count("sim.max_time")) ) {
        std::cerr << "Must provide either sim.max_steps or sim.max_time\n";
        std::exit(1);
    }
}



int main(int argc, const char *argv[])
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
    do {
        double dt = 1e7;
        std::cout << "Step: " << steps << ", time:" << time << "\n";
        steps++;
        time += dt;
    } while (steps <= param.sim.max_steps && time <= param.sim.max_time);


    //
    // output
    //


    return 0;
}
