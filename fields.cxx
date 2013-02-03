#include <iostream>
#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "matprops.hpp"
#include "fields.hpp"


void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    var.edvoldt = new double_vec(e);

    var.temperature = new double_vec(n);
    var.plstrain = new double_vec(e);
    var.ntmp= new double_vec(n);

    var.vel = new array_t(n, 0);
    var.force = new array_t(n, 0);

    var.strain_rate = new tensor_t(e, 0);
    var.spin_rate = new tensor_t(e, 0);
    var.strain = new tensor_t(e, 0);
    var.stress = new tensor_t(e, 0);

    var.shpdx = new shapefn(e);
    if (NDIMS == 3) var.shpdy = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);
}


void reallocate_variables(const Param& param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    delete var.volume;
    delete var.volume_old;
    delete var.volume_n;
    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    delete var.mass;
    delete var.tmass;
    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    delete var.edvoldt;
    var.edvoldt = new double_vec(e);

    delete var.ntmp;
    var.ntmp = new double_vec(n);

    delete var.force;
    var.force = new array_t(n, 0);

    delete var.strain_rate;
    var.strain_rate = new tensor_t(e, 0);

    delete var.shpdx;
    delete var.shpdz;
    var.shpdx = new shapefn(e);
    if (NDIMS == 3) {
        delete var.shpdy;
        var.shpdy = new shapefn(e);
    }
    var.shpdz = new shapefn(e);

    delete var.mat;
    var.mat = new MatProps(param, var);
}


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



void update_spin_rate(const Variables& var, tensor_t& spin_rate)
{
    double *v[NODES_PER_ELEM];

    #pragma omp parallel for default(none) \
        shared(var, spin_rate) private(v)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.spin_rate)[e];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = (*var.vel)[conn[i]];

        // XX component
        int n = 0;
        s[n] = 0;

#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
        // YY component
        n = 1;
        s[n] = 0;
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] - v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] - v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] - v[i][2] * shpdy[i]);
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


void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain)
{

#ifdef THREED
    // The spin rate tensor, W, and the Cauchy stress tensor, S, are         	
    //     [  0  w3  w4]     [s0 s3 s4]
    // W = [-w3   0  w5],  S=[s3 s1 s5].
    //     [-w4 -w5   0]     [s4 s5 s2]
    //
    // Stressi(and strain) increment based on the Jaumann rate is
    // dt*(S*W-W*S). S*W-W*S is also symmetric. 
    // So, following the indexing of stress tensor,
    // we get
    // sj[0] = dt * ( -2 * s3 * w3 - 2 * s4 * w4)
    // sj[1] = dt * (  2 * s3 * w3 - 2 * s5 * w5)
    // sj[2] = dt * (  2 * s4 * w4 + 2 * s5 * w5)
    // sj[3] = dt * ( s0 * w3 - s1 * w3 - s4 * w5 - s5 * w4)
    // sj[4] = dt * ( s0 * w4 - s2 * w4 + s3 * w5 - s5 * w3)
    // sj[5] = dt * ( s1 * w5 - s2 * w5 + s3 * w4 + s4 * w3)

    #pragma omp parallel for default(none) \
        shared(var, stress, strain)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        const double *w = (*var.spin_rate)[e];
        double *s = stress[e];
        double *es = strain[e];
        double s_inc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double es_inc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        s_inc[0] += var.dt * ( -2.0 * s[3] * w[3] - 2.0 * s[4] * w[4] );
        s_inc[1] += var.dt * (  2.0 * s[3] * w[3] - 2.0 * s[5] * w[5] );
        s_inc[2] += var.dt * (  2.0 * s[4] * w[4] + 2.0 * s[5] * w[5] );
        s_inc[3] += var.dt * ( s[0] * w[3] - s[1] * w[3] - s[4] * w[5] - s[5] * w[4] );
        s_inc[4] += var.dt * ( s[0] * w[4] - s[2] * w[4] + s[3] * w[5] - s[5] * w[3] );
        s_inc[5] += var.dt * ( s[1] * w[5] - s[2] * w[5] + s[3] * w[4] + s[4] * w[3] );

        es_inc[0] += var.dt * ( -2.0 * es[3] * w[3] - 2.0 * es[4] * w[4] );
        es_inc[1] += var.dt * (  2.0 * es[3] * w[3] - 2.0 * es[5] * w[5] );
        es_inc[2] += var.dt * (  2.0 * es[4] * w[4] + 2.0 * es[5] * w[5] );
        es_inc[3] += var.dt * ( es[0] * w[3] - es[1] * w[3] - es[4] * w[5] - es[5] * w[4] );
        es_inc[4] += var.dt * ( es[0] * w[4] - es[2] * w[4] + es[3] * w[5] - es[5] * w[3] );
        es_inc[5] += var.dt * ( es[1] * w[5] - es[2] * w[5] + es[3] * w[4] + es[4] * w[3] );

        for(int i=0; i<6; ++i)  {
            s[i] += s_inc[i];
            es[i] += es_inc[i];
        }
    }

#else

    #pragma omp parallel for default(none) \
        shared(var, stress, strain)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double *s = stress[e];
        double *es = strain[e];
        double w = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int nv = conn[i];
            int nx1 = conn[(i+1) % NODES_PER_ELEM];
            int nx2 = conn[(i+2) % NODES_PER_ELEM];
            const double *v = (*var.vel)[nv];
            const double *x1 = (*var.coord)[nx1];
            const double *x2 = (*var.coord)[nx2];
            for (int d=0; d<NDIMS; d++)
                w += v[d] * (x2[d] - x1[d]);
        }

        w *= var.dt * 0.5 / (*var.volume)[e];

        double s12 = s[2];
        double sdiff = s[1] - s[0];
        s[0] += w * s12;
        s[1] -= w * s12;
        s[2] += 0.5 * w * sdiff;

        double e12 = es[2];
        double ediff = es[1] - es[0];
        es[0] += w * e12;
        es[1] -= w * e12;
        es[2] += 0.5 * w * ediff;
    }
#endif
}

