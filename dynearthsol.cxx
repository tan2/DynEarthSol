#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static void allocate_variables(Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    var.dvoldt = new double_vec(n);
    var.edvoldt = new double_vec(e);

    var.temperature = new double_vec(n);
    var.plstrain = new double_vec(e);
    var.tmp0 = new double_vec(std::max(n,e));

    var.vel = new double2d(boost::extents[n][NDIMS]);
    var.force = new double2d(boost::extents[n][NDIMS]);

    var.strain_rate = new double2d(boost::extents[e][NSTR]);
    var.strain = new double2d(boost::extents[e][NSTR]);
    var.stress = new double2d(boost::extents[e][NSTR]);

    var.shpdx = new double2d(boost::extents[e][NODES_PER_ELEM]);
    if (NDIMS == 3) var.shpdy = new double2d(boost::extents[e][NODES_PER_ELEM]);
    var.shpdz = new double2d(boost::extents[e][NODES_PER_ELEM]);
}


static void create_matprops(const Param &par, Variables &var)
{
    // TODO: get material properties from cfg file
    var.mat = new MatProps(par);
}


void compute_dvoldt(const Variables &var, double_vec &tmp,
                    double_vec &dvoldt, double_vec &edvoldt)
{
    /* dvoldt is the volumetric strain rate */
    /* edvoldt is the averaged dvoldt on the element */
    const double_vec& volume = *var.volume;
    const double_vec& volume_n = *var.volume_n;
    std::fill_n(tmp.begin(), var.nnode, 0);

    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                  \
            shared(egroup, var, tmp, volume)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = &(*var.connectivity)[e][0];
        const double* strain_rate = &(*var.strain_rate)[e][0];
        // TODO: try another definition:
        // dj = (volume[e] - volume_old[e]) / volume_old[e] / dt
        double dj = trace(strain_rate);
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            tmp[n] += dj * volume[e];
        }
    }
        } // end of ee
    }


    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, tmp, volume_n)
    for (int n=0; n<var.nnode; ++n)
         dvoldt[n] = tmp[n] / volume_n[n];

    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, edvoldt)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        double dj = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dj += dvoldt[n];
        }
        edvoldt[e] = dj / NODES_PER_ELEM;
    }

    // std::cout << "dvoldt:\n";
    // print(std::cout, dvoldt);
    // std::cout << "\n";
    // std::cout << "edvoldt:\n";
    // print(std::cout, edvoldt);
    // std::cout << "\n";
}


void initial_stress_state(const Param &param, const Variables &var,
                          double2d &stress, double2d &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    // XXX: compute reference pressure correctly
    const double rho = var.mat->density(0);
    const double ks = var.mat->bulkm(0);
    compensation_pressure = rho * param.control.gravity * param.mesh.zlength;
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = rho * param.control.gravity * zcenter;
            strain[e][i] = rho * param.control.gravity * zcenter / ks / NDIMS;
        }
    }
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        // the coordinate of the center of this element
        double center[3] = {0,0,0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }
        // TODO: adding different types of weak zone
        const double d = param.mesh.resolution;
        if (std::fabs(center[0] - param.mesh.xlength * 0.5) < 2*d &&
            std::fabs(center[NDIMS-1] + param.mesh.zlength) < 1.5*d)
            plstrain[e] = 0.1;
    }
}


void initial_temperature(const Param &param, const Variables &var, double_vec &temperature)
{
    const double oceanic_plate_age = 60e6 * YEAR2SEC;
    const double diffusivity = 1e-6;

    for (int i=0; i<var.nnode; ++i) {
        double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * oceanic_plate_age);
        temperature[i] = param.bc.surface_temperature +
            (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
    }
}


void apply_vbcs(const Param &param, const Variables &var, double2d &vel)
{
    // TODO: adding different types of vbcs later

    // diverging x-boundary
    #pragma omp parallel for default(none) \
        shared(param, var, vel)
    for (int i=0; i<var.nnode; ++i) {
        int flag = (*var.bcflag)[i];
        double *v = &vel[i][0];

        // X
        if (flag & BOUNDX0) {
            v[0] = -param.bc.max_vbc_val;
            v[1] = 0;
            v[NDIMS-1] = 0;
        }
        else if (flag & BOUNDX1) {
            v[0] = param.bc.max_vbc_val;
            v[1] = 0;
            v[NDIMS-1] = 0;
        }
#ifdef THREED
        // Y
        if (flag & BOUNDY0) {
            v[1] = 0;
        }
        else if (flag & BOUNDY1) {
            v[1] = 0;
        }
#endif
        // Z
        if (flag & BOUNDZ0) {
            //v[NDIMS-1] = 0;
        }
        else if (flag & BOUNDZ1) {
            //v[NDIMS-1] = 0;
        }
    }
}


void init(const Param& param, Variables& var)
{
    void create_matprops(const Param&, Variables&);

    create_new_mesh(param, var);
    allocate_variables(var);
    create_matprops(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    *var.volume_old = *var.volume;
    compute_dvoldt(var, *var.tmp0, *var.dvoldt, *var.edvoldt);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    initial_stress_state(param, var, *var.stress, *var.strain, var.compensation_pressure);
    initial_weak_zone(param, var, *var.plstrain);
    initial_temperature(param, var, *var.temperature);
    apply_vbcs(param, var, *var.vel);
}


void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{
    // diffusion matrix
    double D[NODES_PER_ELEM][NODES_PER_ELEM];

    tdot.assign(var.nnode, 0);
    for (auto egroup : *var.egroups) {
       #pragma omp parallel for default(none)                                  \
           shared(egroup, var, param, temperature, tdot) private(D)
       for (int ee=0; ee<egroup.size(); ++ee) {
            int e = egroup[ee];
    {
        const int *conn = &(*var.connectivity)[e][0];
        double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volumn
        const double *shpdx = &(*var.shpdx)[e][0];
#ifdef THREED
        const double *shpdy = &(*var.shpdy)[e][0];
#endif
        const double *shpdz = &(*var.shpdz)[e][0];
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


void update_strain_rate(const Variables& var, double2d& strain_rate)
{
    double *v[NODES_PER_ELEM];

    #pragma omp parallel for default(none) \
        shared(var, strain_rate) private(v)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        const double *shpdx = &(*var.shpdx)[e][0];
        const double *shpdz = &(*var.shpdz)[e][0];
        double *s = &(*var.strain_rate)[e][0];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = &(*var.vel)[conn[i]][0];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        const double *shpdy = &(*var.shpdy)[e][0];
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


void update_force(const Param& param, const Variables& var, double2d& force)
{
    std::fill_n(force.data(), var.nnode, 0);

    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                  \
            shared(egroup, param, var, force)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = &(*var.connectivity)[e][0];
        const double *shpdx = &(*var.shpdx)[e][0];
#ifdef THREED
        const double *shpdy = &(*var.shpdy)[e][0];
#endif
        const double *shpdz = &(*var.shpdz)[e][0];
        double *s = &(*var.stress)[e][0];
        double vol = (*var.volume)[e];

        double buoy = var.mat->density(e) * param.control.gravity / NODES_PER_ELEM;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double *f = &force[conn[i]][0];
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


    // Wrinkler foundation
    if (param.bc.wrinkler_foundation && param.control.gravity != 0) {
        const int bottom_bdry = bdry_order.find(BOUNDZ0)->second;
        const auto& bdry = var.bfacets[bottom_bdry];
        const auto& coord = *var.coord;
        for (int i=0; i<bdry.size(); ++i) {
            int e = bdry[i].first;
            const int *conn = &(*var.connectivity)[e][0];
            int facet = bdry[i].second;

            int n0 = conn[NODE_OF_FACET[facet][0]];
            int n1 = conn[NODE_OF_FACET[facet][1]];
#ifdef THREED
            int n2 = conn[NODE_OF_FACET[facet][2]];

            // vectors: n0-n1 and n0-n2
            double v01[NDIMS], v02[NDIMS];
            for (int i=0; i<NDIMS; ++i) {
                v01[i] = coord[n1][i] - coord[n0][i];
                v02[i] = coord[n2][i] - coord[n0][i];
            }

            // the outward normal vector is parallel to the cross product
            // of v01 and v02, with magnitude = area of the triangle
            double normal[NDIMS];
            normal[0] = (v01[1] * v02[2] - v01[2] * v02[1]) / 2;
            normal[1] = (v01[2] * v02[0] - v01[0] * v02[2]) / 2;
            normal[2] = (v01[0] * v02[1] - v01[1] * v02[0]) / 2;

            double zcenter = (coord[n0][2] + coord[n1][2] + coord[n2][2]) / NODES_PER_FACET;
#else
            // the edge vector
            double v01[NDIMS];
            for (int i=0; i<NDIMS; ++i) {
                v01[i] = coord[n1][i] - coord[n0][i];
            }

            // the normal vector to the edge, pointing outward
            double normal[NDIMS];
            normal[0] = v01[1];
            normal[1] = -v01[0];

            double zcenter = (coord[n0][1] + coord[n1][1]) / NODES_PER_FACET;
#endif
            double dz = zcenter - (-param.mesh.zlength);
            double p = var.compensation_pressure -
                (var.mat->density(e) + param.bc.wrinkler_delta_rho) * param.control.gravity * dz;

            // bottom support - Archimed force (normal to the surface)
            for (int i=0; i<NDIMS; ++i) {
                force[n0][i] -= p * normal[i] / NODES_PER_FACET;
                force[n1][i] -= p * normal[i] / NODES_PER_FACET;
#ifdef THREED
                force[n2][i] -= p * normal[i] / NODES_PER_FACET;
#endif
            }
        }
    }


    // damping
    {
        // flatten 2d arrays to simplify indexing
        double* ff = force.data();
        const double* v = var.vel->data();
        const double small_vel = 1e-13;
        const double damping_coeff = 0.8;
        #pragma omp parallel for default(none)  \
            shared(var, ff, v)
        for (int i=0; i<var.nnode*NDIMS; ++i) {
            if (std::fabs(v[i]) > small_vel) {
                ff[i] -= damping_coeff * std::copysign(ff[i], v[i]);
            }
        }
    }
}


void rotate_stress() {};


void update_velocity(const Variables& var, double2d& vel)
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


static void update_coordinate(const Variables& var, double2d_ref& coord)
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


void update_mesh(const Param& param, Variables& var)
{
    update_coordinate(var, *var.coord);

    var.volume->swap(*var.volume_old);
    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    compute_dvoldt(var, *var.tmp0, *var.dvoldt, *var.edvoldt);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);
}


int main(int argc, const char* argv[])
{
    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        std::cout << "       " << argv[0] << " -h or --help\n";
        return -1;
    }

    Param param;
    void get_input_parameters(const char*, Param&);
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    var.time = 0;
    var.steps = 0;
    var.frame = 0;

    if (! param.sim.is_restarting) {
        init(param, var);
        output(param, var);
        var.frame ++;
    }
    else {
        restart();
        var.frame ++;
    }

    var.dt = compute_dt(param, var);

    do {
        var.steps ++;
        var.time += var.dt;

        update_temperature(param, var, *var.temperature, *var.tmp0);
        update_strain_rate(var, *var.strain_rate);
        update_stress(var, *var.stress, *var.strain, *var.plstrain, *var.strain_rate);
        update_force(param, var, *var.force);
        update_velocity(var, *var.vel);
        apply_vbcs(param, var, *var.vel);
        update_mesh(param, var);
        // dt computation is expensive, and dt only changes slowly
        // don't have to do it every time step
        if (var.steps % 10 == 0) var.dt = compute_dt(param, var);
        rotate_stress();

        if ( (var.steps == var.frame * param.sim.output_step_interval) ||
             (var.time > var.frame * param.sim.output_time_interval_in_yr * YEAR2SEC) ) {
            output(param, var);
            var.frame ++;
        }

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    return 0;
}
