
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

    po::options_description cfg("Config file options");
    cfg.add_options()
        ("sim.steps", po::value<int>(), "Max. number of time steps")
        ;

    cfg.add_options()
        ("mesh.xlength", po::value<double>(), "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(), "Length of y (in meters)")
        ("mesh.zlength", po::value<double>(), "Length of z (in meters)")
        ;

    //
    // read config file
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
    // get parameters
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


    //
    // output
    //


    return 0;
}
