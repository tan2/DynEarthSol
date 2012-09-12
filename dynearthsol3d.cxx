
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
//#include <string>

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
    po::options_description cfg("Config file options");
    cfg.add_options()
        ("sim.max_steps", po::value<int>(), "Max. number of time steps")
        ("sim.max_time", po::value<double>(), "Max. time (in seconds)")
        ;

    cfg.add_options()
        ("mesh.xlength", po::value<double>(), "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(), "Length of y (in meters)")
        ("mesh.zlength", po::value<double>(), "Length of z (in meters)")
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
    // get parameters from config file
    //
    try {
        if (vm.count("mesh.xlength")) {
            std::cout << vm["mesh.xlength"].as<double>() << "\n";
        }
    }
    catch (std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    //
    // run simulation
    //
    int steps = 0;
    double time = 0;
    int max_steps = vm["sim.max_steps"].as<int>();
    double max_time = vm["sim.max_time"].as<double>();
    do {
        double dt = 0;
        std::cout << "Step: " << steps << ", time:" << max_time << "\n";
        steps++;
        time += dt;
    } while (steps <= max_steps && time <= max_time);


    //
    // output
    //


    return 0;
}
