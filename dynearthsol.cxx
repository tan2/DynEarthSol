#include <cstdio>
#include <ctime>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
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

    var.jacobian = new double_vec(n);
    var.ejacobian = new double_vec(e);

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
    var.mat = new MatProps(1, MatProps::rh_evp);
}


void initial_stress_state(const Param &param, const Variables &var,
                          double2d &stress, double2d &strain,
                          double &compensation_pressure)
{
    if (param.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    // XXX: compute reference pressure correctly
    const double rho = var.mat->density(0);
    const double ks = var.mat->bulkm(0);
    compensation_pressure = rho * param.gravity * param.mesh.zlength;
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = - rho * param.gravity * zcenter;
            strain[e][i] = - rho * param.gravity * zcenter / ks / NDIMS;
        }
    }
}


void initial_temperature(const Param &param, const Variables &var, double_vec &temperature)
{
    const double oceanic_plate_age = 1e6 * YEAR2SEC;
    const double diffusivity = 1e-6;

    for (int i=0; i<var.nnode; ++i) {
        double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * oceanic_plate_age);
        temperature[i] = param.surface_temperature +
            (param.mantle_temperature - param.surface_temperature) * std::erf(w);
    }
}


void apply_vbcs(const Param &param, const Variables &var, double2d &vel)
{
    // TODO: adding different types of vbcs later

    // diverging x-boundary
    for (int i=0; i<var.nnode; ++i) {
        int flag = (*var.bcflag)[i];

        // X
        if (flag & BOUNDX0) {
            vel[i][0] = -param.maxvbcval;
        }
        else if (flag & BOUNDX1) {
            vel[i][0] = param.maxvbcval;
        }
#ifdef THREED
        // Y
        if (flag & BOUNDY0) {
            vel[i][1] = 0;
        }
        else if (flag & BOUNDY1) {
            vel[i][1] = 0;
        }
#endif
        // Z
        if (flag & BOUNDZ0) {
            //vel[i][NDIMS-1] = 0;
        }
        else if (flag & BOUNDZ1) {
            vel[i][NDIMS-1] = 0;
        }
    }
}


void init(const Param& param, Variables& var)
{
    void create_matprops(const Param&, Variables&);

    create_new_mesh(param, var);
    allocate_variables(var);
    create_matprops(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.volume, *var.volume_n);
    compute_mass(*var.coord, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume,
                     *var.shpdx, *var.shpdy, *var.shpdz);
    // XXX
    //create_jacobian();

    initial_stress_state(param, var, *var.stress, *var.strain, var.compensation_pressure);
    initial_temperature(param, var, *var.temperature);
    apply_vbcs(param, var, *var.vel);
};


void restart() {};


void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{
    // diffusion matrix
    double D[NODES_PER_ELEM][NODES_PER_ELEM];

    tdot.assign(var.nnode, 0);
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = &(*var.connectivity)[e][0];
        double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volumn
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int j=0; j<NODES_PER_ELEM; ++j) {
                if (NDIMS == 3) {
                    D[i][j] = ((*var.shpdx)[e][i] * (*var.shpdx)[e][j] +
                               (*var.shpdy)[e][i] * (*var.shpdy)[e][j] +
                               (*var.shpdz)[e][i] * (*var.shpdz)[e][j]);
                }
                else {
                    D[i][j] = ((*var.shpdx)[e][i] * (*var.shpdx)[e][j] +
                               (*var.shpdz)[e][i] * (*var.shpdz)[e][j]);
                }
            }
        }
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double diffusion = 0;
            for (int j=0; j<NODES_PER_ELEM; ++j)
                diffusion += D[i][j] * temperature[conn[j]];

            tdot[conn[i]] += diffusion * kv;
        }
    }

    for (int n=0; n<var.nnode; ++n) {
        if ((*var.bcflag)[n] & BOUNDZ1)
            temperature[n] = param.surface_temperature;
        else
            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
    }
}


void update_strain_rate() {};
void update_stress() {};
void update_force() {};
void rotate_stress() {};


static void update_coordinate(const Variables& var, double2d_ref& coord)
{
    double* x = var.coord->data();
    const double* v = var.vel->data();
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        x[i] += v[i] * var.dt;
    }

    // surface_processes()
}


void update_mesh(const Param& param, Variables& var)
{
    update_coordinate(var, *var.coord);

    var.volume->swap(*var.volume_old);
    compute_volume(*var.coord, *var.connectivity, *var.volume, *var.volume_n);
    compute_mass(*var.coord, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume,
                     *var.shpdx, *var.shpdy, *var.shpdz);
};


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

    // temperature
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "temperature", var.frame);
    write_array(buffer, *var.temperature);
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
        update_strain_rate();
        update_stress();
        update_force();
        update_mesh(param, var);
        // dt computation is expensive, and dt only changes slowly
        // don't have to do it every time step
        if (var.steps % 10 == 0) var.dt = compute_dt(param, var);
        rotate_stress();

        if ( (var.steps == var.frame * param.sim.output_step_interval) ||
             (var.time > var.frame * param.sim.output_time_interval_in_yr * YEAR2SEC) ) {
            output(param, var);
            std::cout << "  Output # " << var.frame
                      << ", step = " << var.steps
                      << ", time = " << var.time / YEAR2SEC << " yr"
                      << ", dt = " << var.dt / YEAR2SEC << " yr.\n";
            var.frame ++;
        }

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    return 0;
}
