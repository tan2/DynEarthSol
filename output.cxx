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

/*****************************************************************************
 * The format of the .save file:
 * 1  The first 'headerlen' bytes are ASCII text.
 *   1.1  The 1st line in the header is the revision string. Starting with
 *        "# DynEarthSol ndims=%1 revision=%2", with %1 equal to 2 or 3
 *        (indicating 2D or 3D simulation) and %2 an integer.
 *   1.2  The following lines are 'name', 'position' pairs, separated by a
 *        TAB character. This line tells the name of the data and the
 *        starting position (in bytes) of the data in this file.
 * 2  The rests are binary data.
 ****************************************************************************/


Output::Output(const Param& param, double start_time, int start_frame) :
    modelname(param.sim.modelname),
    start_time(start_time),
    average_interval(param.sim.output_averaged_fields),
    frame(start_frame),
    time0(0)
{}


Output::~Output()
{}


template <typename T>
void Output::write_array(std::FILE* f, const std::vector<T>& A, const char *name)
{
    write_header(f, name);
    std::size_t n = std::fwrite(A.data(), sizeof(T), A.size(), f);
    eof_pos += n * sizeof(T);
}


template <typename T, int N>
void Output::write_array(std::FILE* f, const Array2D<T,N>& A, const char *name)
{
    write_header(f, name);
    std::size_t n = std::fwrite(A.data(), sizeof(T), A.num_elements(), f);
    eof_pos += n * sizeof(T);
}


void Output::write_header(std::FILE* f, const char *name)
{
    const int bsize = 256;
    char buffer[bsize];
    int len = std::snprintf(buffer, bsize, "%s\t%ld\n", name, eof_pos);
    if (len >= bsize) {
        std::cerr << "Error: exceeding buffer length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        std::exit(1);
    }
    if (len >= headerlen - (hd_pos - header)*sizeof(char)) {
        std::cerr << "Error: exceeding header length at Output::write_array, name=" << name
                  << " eof_position=" << eof_pos << '\n';
        std::exit(1);
    }
    hd_pos = std::strncat(hd_pos, buffer, len);
}


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
    char filename[256];
    std::FILE* f;

    double dt = var.dt;
    if (average_interval && is_averaged)
        dt = (var.time - time0) / average_interval;

    write_info(var, dt);

    std::snprintf(filename, 255, "%s.save.%06d", modelname.c_str(), frame);
    f = fopen(filename, "w");

    header = new char[headerlen]();
    hd_pos = std::strcat(header, revision_str);
    eof_pos = headerlen;
    std::fseek(f, eof_pos, SEEK_SET);

    write_array(f, *var.coord, "coordinate");
    write_array(f, *var.connectivity, "connectivity");

    array_t *vel = var.vel;
    if (average_interval && is_averaged) {
        // average_velocity = displacement / delta_t
        vel = &coord0;
        double *c0 = coord0.data();
        const double *c = var.coord->data();
        for (int i=0; i<coord0.num_elements(); ++i) {
            c0[i] = (c[i] - c0[i]) / (var.time - time0);
        }
    }
    write_array(f, *vel, "velocity");

    write_array(f, *var.temperature, "temperature");


    write_array(f, *var.elquality, "mesh quality");
    write_array(f, *var.plstrain, "plastic strain");

    double_vec *delta_plstrain = var.delta_plstrain;
    if (average_interval && is_averaged) {
        delta_plstrain = &delta_plstrain_avg;
        double tmp = 1.0 / (average_interval + 1);
        for (int i=0; i<delta_plstrain_avg.size(); ++i) {
            delta_plstrain_avg[i] *= tmp;
        }
    }
    write_array(f, *delta_plstrain, "plastic strain increment");

    tensor_t *strain_rate = var.strain_rate;
    if (average_interval && is_averaged) {
        // average_strain_rate = delta_strain / delta_t
        strain_rate = &strain0;
        double *s0 = strain0.data();
        const double *s = var.strain->data();
        for (int i=0; i<strain0.num_elements(); ++i) {
            s0[i] = (s[i] - s0[i]) / (var.time - time0);
        }
    }
    write_array(f, *strain_rate, "strain-rate");

    write_array(f, *var.strain, "strain");

    tensor_t *stress = var.stress;
    if (average_interval && is_averaged) {
        stress = &stress_avg;
        double *s = stress_avg.data();
        double tmp = 1.0 / (average_interval + 1);
        for (int i=0; i<stress_avg.num_elements(); ++i) {
            s[i] *= tmp;
        }
    }
    write_array(f, *stress, "stress");

    double_vec tmp(var.nelem);
    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->rho(e);
    }
    write_array(f, tmp, "density");

    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->visc(e);
    }
    write_array(f, tmp, "viscosity");

    write_array(f, *var.volume, "volume");
    write_array(f, *var.volume_old, "volume_old");

    write_array(f, *var.force, "force");

    //write_array(f, *var.bcflag, "bcflag");

    /* write header */
    std::fseek(f, 0, SEEK_SET);
    std::fwrite(header, sizeof(char), headerlen, f);
    delete [] header;

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


