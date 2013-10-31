#include <cstdio>
#include <iostream>

#ifdef USE_OMP
#include <omp.h>
#else
#include <ctime>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"


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


void output(const Param& param, const Variables& var, double start_time)
{
    /* Not using C++ stream IO here since it can be much slower than C stdio. */

    using namespace std;
    char buffer[255];
    std::FILE* f;

#ifdef USE_OMP
    double run_time = omp_get_wtime() - start_time;
#else
    double run_time = double(std::clock()) / CLOCKS_PER_SEC;
#endif

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

    // mesh quality
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "meshquality", var.frame);
    write_array(buffer, *var.elquality);

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

    // density
    double_vec tmp(var.nelem);
    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->rho(e);
    }
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "density", var.frame);
    write_array(buffer, tmp);

    // viscosity
    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->visc(e);
    }
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "viscosity", var.frame);
    write_array(buffer, tmp);

    // volume
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "volume", var.frame);
    write_array(buffer, *var.volume);

    // volume_old
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "volume_old", var.frame);
    write_array(buffer, *var.volume_old);

    // force with boundary removed
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "force", var.frame);
    write_array(buffer, *var.force);

    // bcflag
    //snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "bcflag", var.frame);
    //write_array(buffer, *var.bcflag);

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


