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
#include "output.hpp"


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


Output::Output(const Param& param, double start_time, int start_frame) :
    modelname(param.sim.modelname),
    start_time(start_time),
    average_interval(param.sim.output_averaged_fields),
    frame(start_frame),
    time0(0)
{}


Output::~Output()
{}


void Output::write_info(const Variables& var, double dt)
{
#ifdef USE_OMP
    double run_time = omp_get_wtime() - start_time;
#else
    double run_time = double(std::clock()) / CLOCKS_PER_SEC;
#endif

    char buffer[256];
    std::snprintf(buffer, 255, "%6d\t%10d\t%12.6e\t%12.4e\t%12.6e\t%8d\t%8d\t%8d\n",
                  frame, var.steps, var.time, dt, run_time,
                  var.nnode, var.nelem, var.nseg);

    std::string filename(modelname + ".info");
    std::FILE* f;
    if (frame == 0)
        f = std::fopen(filename.c_str(), "w");
    else
        f = std::fopen(filename.c_str(), "a");
    std::fputs(buffer, f);
    std::fclose(f);
}


void Output::write(const Variables& var, bool is_averaged)
{
    /* Not using C++ stream IO here since it can be much slower than C stdio. */

    using namespace std;
    char buffer[255];
    std::FILE* f;

    double dt = var.dt;
    if (average_interval && is_averaged)
        dt = (var.time - time0) / average_interval;

    write_info(var, dt);

    // coord
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "coord", frame);
    write_array(buffer, *var.coord);

    // connectivity
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "connectivity", frame);
    write_array(buffer, *var.connectivity);

    // vel
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "vel", frame);
    if (average_interval && is_averaged) {
        // average_velocity = displacement / delta_t
        double *c0 = coord0.data();
        const double *c = var.coord->data();
        for (int i=0; i<coord0.num_elements(); ++i) {
            c0[i] = (c[i] - c0[i]) / (var.time - time0);
        }
        write_array(buffer, coord0);
    }
    else
        write_array(buffer, *var.vel);

    // temperature
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "temperature", frame);
    write_array(buffer, *var.temperature);

    // mesh quality
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "meshquality", frame);
    write_array(buffer, *var.elquality);

    // plstrain
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "plstrain", frame);
    write_array(buffer, *var.plstrain);

    // delta_plstrain
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "delta_plstrain", frame);
    if (average_interval && is_averaged) {
        double tmp = 1.0 / (average_interval + 1);
        for (int i=0; i<delta_plstrain_avg.size(); ++i) {
            delta_plstrain_avg[i] *= tmp;
        }
        write_array(buffer, delta_plstrain_avg);
    }
    else
        write_array(buffer, *var.delta_plstrain);

    // strain_rate
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "strain-rate", frame);
    if (average_interval && is_averaged) {
        // average_strain_rate = delta_strain / delta_t
        double *s0 = strain0.data();
        const double *s = var.strain->data();
        for (int i=0; i<strain0.num_elements(); ++i) {
            s0[i] = (s[i] - s0[i]) / (var.time - time0);
        }
        write_array(buffer, strain0);
    }
    else
        write_array(buffer, *var.strain_rate);

    // strain
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "strain", frame);
    write_array(buffer, *var.strain);

    // stress
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "stress", frame);
    if (average_interval && is_averaged) {
        double *s = stress_avg.data();
        double tmp = 1.0 / (average_interval + 1);
        for (int i=0; i<stress_avg.num_elements(); ++i) {
            s[i] *= tmp;
        }
        write_array(buffer, stress_avg);
    }
    else
        write_array(buffer, *var.stress);

    // density
    double_vec tmp(var.nelem);
    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->rho(e);
    }
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "density", frame);
    write_array(buffer, tmp);

    // viscosity
    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->visc(e);
    }
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "viscosity", frame);
    write_array(buffer, tmp);

    // volume
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "volume", frame);
    write_array(buffer, *var.volume);

    // volume_old
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "volume_old", frame);
    write_array(buffer, *var.volume_old);

    // force with boundary removed
    snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "force", frame);
    write_array(buffer, *var.force);

    // bcflag
    //snprintf(buffer, 255, "%s.%s.%06d", modelname.c_str(), "bcflag", frame);
    //write_array(buffer, *var.bcflag);

    std::cout << "  Output # " << frame
              << ", step = " << var.steps
              << ", time = " << var.time / YEAR2SEC << " yr"
              << ", dt = " << dt / YEAR2SEC << " yr.\n";

    frame ++;
}


void Output::average_fields(Variables& var)
{
    // In the first time step of each interval, some old fields are
    // stored for later use.
    // It is guaranteed that remeshing won't occur within the interval.
    if (var.steps % average_interval == 1) {
        time0 = var.time;

        if (coord0.size() != var.coord->size()) {
            double *tmp = new double[(var.coord)->num_elements()];
            coord0.reset(tmp, var.coord->size());
        }
        std::copy(var.coord->begin(), var.coord->end(), coord0.begin());

        if (strain0.size() != var.strain->size()) {
            double *tmp = new double[(var.strain)->num_elements()];
            strain0.reset(tmp, var.strain->size());
        }
        std::copy(var.strain->begin(), var.strain->end(), strain0.begin());

        if (stress_avg.size() != var.stress->size()) {
            double *tmp = new double[(var.stress)->num_elements()];
            stress_avg.reset(tmp, var.stress->size());
        }
        std::copy(var.stress->begin(), var.stress->end(), stress_avg.begin());

        delta_plstrain_avg = *var.delta_plstrain;
    }
    else {
        // Averaging stress & plastic strain
        // (PS: dt-weighted average would be better, but difficult to do)
        double *s_avg = stress_avg.data();
        const double *s = var.stress->data();
        for (int i=0; i<stress_avg.num_elements(); ++i) {
            s_avg[i] += s[i];
        }
        for (int i=0; i<delta_plstrain_avg.size(); ++i) {
            delta_plstrain_avg[i] += (*var.delta_plstrain)[i];
        }
    }
}


void restart()
{
    std::cerr << "Error: restarting not implemented.\n";
    std::exit(1);

    std::cout << "Initializing mesh and field data from checkpoints...\n";
}


