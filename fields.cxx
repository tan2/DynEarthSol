#include <iostream>
#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "matprops.hpp"



void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{
    // diffusion matrix
    double D[NODES_PER_ELEM][NODES_PER_ELEM];

    tdot.assign(var.nnode, 0);
    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                                  \
            shared(egroup, var, param, temperature, tdot) private(D)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = (*var.connectivity)[e];
        double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volumn
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int j=0; j<NODES_PER_ELEM; ++j) {
#ifdef THREED
                D[i][j] = (shpdx[i] * shpdx[j] +
                           shpdy[i] * shpdy[j] +
                           shpdz[i] * shpdz[j]);
#else
                D[i][j] = (shpdx[i] * shpdx[j] +
                           shpdz[i] * shpdz[j]);
#endif
            }
        }
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double diffusion = 0;
            for (int j=0; j<NODES_PER_ELEM; ++j)
                diffusion += D[i][j] * temperature[conn[j]];

            tdot[conn[i]] += diffusion * kv;
        }
    }
        }
    }

     #pragma omp parallel for default(none)      \
         shared(var, param, tdot, temperature)
     for (int n=0; n<var.nnode; ++n) {
        if ((*var.bcflag)[n] & BOUNDZ1)
            temperature[n] = param.bc.surface_temperature;
        else
            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
    }
}


void update_strain_rate(const Variables& var, tensor_t& strain_rate)
{
    double *v[NODES_PER_ELEM];

    #pragma omp parallel for default(none) \
        shared(var, strain_rate) private(v)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.strain_rate)[e];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = (*var.vel)[conn[i]];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
        // YY component
        n = 1;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][1] * shpdy[i];
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][NDIMS-1] * shpdz[i];

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] + v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] + v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] + v[i][2] * shpdy[i]);
#endif
    }
}


static void apply_damping(const Param& param, const Variables& var, array_t& force)
{
    // flatten 2d arrays to simplify indexing
    double* ff = force.data();
    const double* v = var.vel->data();
    const double small_vel = 1e-13;
    #pragma omp parallel for default(none)          \
        shared(var, param, ff, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        if (std::fabs(v[i]) > small_vel) {
            ff[i] -= param.control.damping_factor * std::copysign(ff[i], v[i]);
        }
    }
}


void update_force(const Param& param, const Variables& var, array_t& force)
{
    std::fill_n(force.data(), var.nnode*NDIMS, 0);

    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                  \
            shared(egroup, param, var, force)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.stress)[e];
        double vol = (*var.volume)[e];

        double buoy = 0;
        if (param.control.gravity != 0)
            buoy = var.mat->rho(e) * param.control.gravity / NODES_PER_ELEM;

        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double *f = force[conn[i]];
#ifdef THREED
            f[0] -= (s[0]*shpdx[i] + s[3]*shpdy[i] + s[4]*shpdz[i]) * vol;
            f[1] -= (s[3]*shpdx[i] + s[1]*shpdy[i] + s[5]*shpdz[i]) * vol;
            f[2] -= (s[4]*shpdx[i] + s[5]*shpdy[i] + s[2]*shpdz[i] + buoy) * vol;
#else
            f[0] -= (s[0]*shpdx[i] + s[2]*shpdz[i]) * vol;
            f[1] -= (s[2]*shpdx[i] + s[1]*shpdz[i] + buoy) * vol;
#endif
        }
    }
        } // end of ee
    }

    apply_stress_bcs(param, var, force);

    if (param.control.damping_factor != 0) {
        apply_damping(param, var, force);
    }
}


void update_velocity(const Variables& var, array_t& vel)
{
    const double* m = &(*var.mass)[0];
    // flatten 2d arrays to simplify indexing
    const double* f = var.force->data();
    double* v = vel.data();
    #pragma omp parallel for default(none) \
        shared(var, m, f, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        int n = i / NDIMS;
        v[i] += var.dt * f[i] / m[n];
    }
}


void update_coordinate(const Variables& var, array_t& coord)
{
    double* x = var.coord->data();
    const double* v = var.vel->data();

    #pragma omp parallel for default(none) \
        shared(var, x, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        x[i] += v[i] * var.dt;
    }

    // surface_processes()
}


void rotate_stress() {};


