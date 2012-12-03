#include <cstdio>
#include <ctime>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"


template <typename Array>
static void write_array(const char* filename, const Array& A)
{
    std::FILE *f = fopen(filename, "w");
    std::fwrite(A.data(), sizeof(typename Array::element), A.num_elements(), f);
    std::fclose(f);
}


template <typename T>
static void write_array(const char* filename, const std::vector<T>& A)
{
    std::FILE *f = fopen(filename, "w");
    std::fwrite(&A.front(), sizeof(T), A.size(), f);
    std::fclose(f);
}


void output(const Param& param, const Variables& var)
{
    /* Not using C++ stream IO here since it can be much slower than C stdio. */

    using namespace std;
    char buffer[255];
    std::FILE* f;

    double run_time = double(std::clock()) / CLOCKS_PER_SEC;

    // info
    snprintf(buffer, 255, "%s.%s", param.sim.modelname.c_str(), "info");
    if (var.frame == 0)
        f = fopen(buffer, "w");
    else
        f = fopen(buffer, "a");

    snprintf(buffer, 255, "%6d\t%10d\t%12.6e\t%12.4e\t%12.6e\t%8d\t%8d\t%8d\n",
             var.frame, var.steps, var.time, var.dt, run_time,
             var.nnode, var.nelem, var.nseg);
    fputs(buffer, f);
    fclose(f);

    // coord
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "coord", var.frame);
    write_array(buffer, *var.coord);

    // connectivity
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "connectivity", var.frame);
    write_array(buffer, *var.connectivity);

    // vel
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "vel", var.frame);
    write_array(buffer, *var.vel);

    // temperature
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "temperature", var.frame);
    write_array(buffer, *var.temperature);

    // plstrain
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "plstrain", var.frame);
    write_array(buffer, *var.plstrain);

    // strain_rate
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "strain-rate", var.frame);
    write_array(buffer, *var.strain_rate);

    // strain
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "strain", var.frame);
    write_array(buffer, *var.strain);

    // stress
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "stress", var.frame);
    write_array(buffer, *var.stress);

    // volume
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "volume", var.frame);
    write_array(buffer, *var.volume);

    // volume_old
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "volume_old", var.frame);
    write_array(buffer, *var.volume_old);

    // force with boundary removed
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "force", var.frame);
    write_array(buffer, *var.force);

    std::cout << "  Output # " << var.frame
              << ", step = " << var.steps
              << ", time = " << var.time / YEAR2SEC << " yr"
              << ", dt = " << var.dt / YEAR2SEC << " yr.\n";
}


void restart()
{
    std::cerr << "Error: restarting not implemented.\n";
    std::exit(1);
}


