#include <cstdio>
#include <limits>
#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "parameters.hpp"
#include "mesh.hpp"

void get_input_parameters(const char* filename, Param& p)
{
    //
    // declare input parameters
    //
    po::options_description cfg("Config file options");
    cfg.add_options()
        ("sim.modelname", po::value<std::string>(&p.sim.modelname), "Prefix for the output files")

        ("sim.max_steps", po::value<int>(&p.sim.max_steps), "Max. number of time steps")
        ("sim.max_time", po::value<double>(&p.sim.max_time), "Max. time (in seconds)")
        ("sim.output_step_interval", po::value<int>(&p.sim.output_step_interval),
         "Output step interval")
        ("sim.output_time_interval", po::value<double>(&p.sim.output_time_interval),
         "Output time interval")
        ("sim.is_restarting", po::value<bool>(&p.sim.is_restarting)->default_value(false),
         "Restarting from previous save?")
        ;

    cfg.add_options()
        ("mesh.xlength", po::value<double>(&p.mesh.xlength)->required(), "Length of x (in meters)")
        ("mesh.ylength", po::value<double>(&p.mesh.ylength)->required(), "Length of y (in meters)")
        ("mesh.zlength", po::value<double>(&p.mesh.zlength)->required(), "Length of z (in meters)")
        ("mesh.resolution", po::value<double>(&p.mesh.resolution)->required(), "Spatial resolution (in meters)")
        ;

    //
    // parse config file
    //
    po::variables_map vm;
    try {
        po::store(po::parse_config_file<char>(filename, cfg), vm);
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
        p.sim.max_steps = std::numeric_limits<int>::max();;
    if ( ! vm.count("sim.max_time") )
        p.sim.max_time = std::numeric_limits<double>::max();;

    if ( ! (vm.count("sim.output_step_interval") || vm.count("sim.output_time_interval")) ) {
        std::cerr << "Must provide either sim.output_step_interval or sim.output_time_interval\n";
        std::exit(1);
    }
    if ( ! vm.count("sim.output_step_interval") )
        p.sim.output_step_interval = std::numeric_limits<int>::max();;
    if ( ! vm.count("sim.output_time_interval") )
        p.sim.output_time_interval = std::numeric_limits<double>::max();;

    return;
}


void init(const Param& param, Variables& var)
{
    new_mesh(param, var);
};


void restart() {};
void update_temperature() {};
void update_strain_rate() {};
void update_stress() {};
void update_force() {};
void update_mesh() {};
void rotate_stress() {};


void output(const Param& param, const Variables& var,
            int steps, int frame, double time, double dt)
{
    /* Not using C++ stream IO here since it can be much slower than C stdio. */

    using namespace std;
    char buffer[255];
    std::FILE* f;

    double run_time = 3.548641321351658e-13; // XXX

    // info
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "info", frame);
    if (frame == 0)
        f = fopen(buffer, "w");
    else
        f = fopen(buffer, "a");

    snprintf(buffer, 255, "%6d\t%10d\t%12.6e\t%12.4e\t%12.6e\t%8d\t%8d\t%8d\n",
                  frame, steps, time, dt, run_time, 0, 0, 0);
    fputs(buffer, f);
    fclose(f);

    // coord
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "coord", frame);
    f = fopen(buffer, "w");
    fwrite(var.coord->data(), sizeof(double), var.coord->num_elements(), f);
    fclose(f);

    // connectivity
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "connectivity", frame);
    f = fopen(buffer, "w");
    fwrite(var.connectivity->data(), sizeof(int), var.connectivity->num_elements(), f);
    fclose(f);

}


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
    Variables var;
    int steps = 0;
    int frame = 0;
    double time = 0;
    double dt = 1e7;
    if (! param.sim.is_restarting) {
        init(param, var);
        output(param, var, steps, frame, time, dt);
        frame ++;
    }
    else {
        restart();
        frame ++;
    }

    do {
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
            //output();
            std::cout << frame <<"-th output\n";
            frame++;
        }

    } while (steps < param.sim.max_steps && time <= param.sim.max_time);

    return 0;
}
