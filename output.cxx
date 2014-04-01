#include <cmath>
#include <cstdio>
#include <iostream>

#ifdef USE_OMP
#include <omp.h>
#else
#include <ctime>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"
#include "markerset.hpp"
#include "matprops.hpp"
#include "output.hpp"


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
    double dt = var.dt;
    double inv_dt = 1 / var.dt;
    if (average_interval && is_averaged) {
        dt = (var.time - time0) / average_interval;
        inv_dt = 1.0 / (var.time - time0);
    }
    write_info(var, dt);

    char filename[256];
    std::snprintf(filename, 255, "%s.save.%06d", modelname.c_str(), frame);
    BinaryOutput bin(filename);

    bin.write_array(*var.coord, "coordinate");
    bin.write_array(*var.connectivity, "connectivity");

    bin.write_array(*var.vel, "velocity");
    if (average_interval && is_averaged) {
        // average_velocity = displacement / delta_t
        double *c0 = coord0.data();
        const double *c = var.coord->data();
        for (int i=0; i<coord0.num_elements(); ++i) {
            c0[i] = (c[i] - c0[i]) * inv_dt;
        }
        bin.write_array(coord0, "velocity averaged");
    }

    bin.write_array(*var.temperature, "temperature");


    bin.write_array(*var.elquality, "mesh quality");
    bin.write_array(*var.plstrain, "plastic strain");

    // Strain rate and plastic strain rate do not need to be checkpointed,
    // so we don't have to distinguish averged/non-averaged variants.
    double_vec *delta_plstrain = var.delta_plstrain;
    if (average_interval && is_averaged) {
        // average_strain_rate = delta_strain / delta_t
        delta_plstrain = &delta_plstrain_avg;
    }
    for (int i=0; i<delta_plstrain_avg.size(); ++i) {
        delta_plstrain_avg[i] *= inv_dt;
    }
    bin.write_array(*delta_plstrain, "plastic strain-rate");

    tensor_t *strain_rate = var.strain_rate;
    if (average_interval && is_averaged) {
        // average_strain_rate = delta_strain / delta_t
        strain_rate = &strain0;
        double *s0 = strain0.data();
        const double *s = var.strain->data();
        for (int i=0; i<strain0.num_elements(); ++i) {
            s0[i] = (s[i] - s0[i]) * inv_dt;
        }
    }
    bin.write_array(*strain_rate, "strain-rate");

    bin.write_array(*var.strain, "strain");
    bin.write_array(*var.stress, "stress");

    if (average_interval && is_averaged) {
        double *s = stress_avg.data();
        double tmp = 1.0 / (average_interval + 1);
        for (int i=0; i<stress_avg.num_elements(); ++i) {
            s[i] *= tmp;
        }
        bin.write_array(stress_avg, "stress averaged");
    }

    double_vec tmp(var.nelem);
    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->rho(e);
    }
    bin.write_array(tmp, "density");

    for (int e=0; e<var.nelem; ++e) {
        tmp[e] = var.mat->visc(e);
    }
    bin.write_array(tmp, "viscosity");

    //bin.write_array(*var.volume, "volume");

    bin.write_array(*var.force, "force");

    //bin.write_array(*var.bcflag, "bcflag");

    bin.close();
    std::cout << "  Output # " << frame
              << ", step = " << var.steps
              << ", time = " << var.time / YEAR2SEC << " yr"
              << ", dt = " << dt / YEAR2SEC << " yr.\n";

    frame ++;

    {
        // check for NaN in coordinate
        for (int i=0; i<var.nnode; i++)
            for (int j=0; j<NDIMS; j++) {
                if (std::isnan((*var.coord)[i][j])) {
                    std::cerr << "Error: coordinate becomes NaN\n";
                    std::exit(1);
                }
                if (std::isinf((*var.coord)[i][j])) {
                    std::cerr << "Error: coordinate becomes Infinity\n";
                    std::exit(1);
                }
            }
    }
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


void Output::write_checkpoint(const Variables& var)
{
    char filename[256];
    std::snprintf(filename, 255, "%s.chkpt.%06d", modelname.c_str(), frame);
    BinaryOutput bin(filename);

    double_vec tmp(2);
    tmp[0] = var.time;
    tmp[1] = var.compensation_pressure;
    bin.write_array(tmp, "time compensation_pressure");

    bin.write_array(*var.segment, "segment");
    bin.write_array(*var.segflag, "segflag");
    // Note: regattr is not needed for restarting
    // bin.write_array(*var.regattr, "regattr");

    bin.write_array(*var.volume_old, "volume_old");

    MarkerSet &ms = *var.markerset;
    ms.write(bin);
}

