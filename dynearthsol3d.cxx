
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
//#include <string>


struct param {
    int max_steps;
    double max_time;

    double xlength, ylength, zlength;
};


int main(int argc, const char *argv[])
{
    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        return -1;
    }

    std::ifstream ifs(argv[1]);
    if (!ifs) {
        std::cerr << "Error: cannot open config_file '" << argv[1] << "'\n";
        return 1;
    }

    //
    // declare input parameters
    //
    struct param param;
    po::options_description cfg("Config file options");
    cfg.add_options()
        ("sim.max_steps", po::value<int>(&param.max_steps), "Max. number of time steps")
        ("sim.max_time", po::value<double>(&param.max_time), "Max. time (in seconds)")
        ;

    cfg.add_options()
        ("mesh.xlength", po::value<double>(&param.xlength)->required(), "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(&param.ylength)->required(), "Length of y (in meters)")
        ("mesh.zlength", po::value<double>(&param.zlength)->required(), "Length of z (in meters)")
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
        std::cerr << "Error reading config_file '" << argv[1] << "'\n";
        std::cerr << e.what() << "\n";
        return 1;
    }

    //
    // validate parameters
    //
    if ( !(vm.count("sim.max_steps") || vm.count("sim.max_time")) ) {
        std::cerr << "Must provide either sim.max_steps or sim.max_time\n";
        return 1;
    }

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
    } while (steps <= param.max_steps && time <= param.max_time);


    //
    // output
    //


    return 0;
}
