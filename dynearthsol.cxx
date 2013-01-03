#include <iostream>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "ic.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "rheology.hpp"


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

    var.vel = new array_t(n, 0);
    var.force = new array_t(n, 0);

    var.strain_rate = new tensor_t(e, 0);
    var.strain = new tensor_t(e, 0);
    var.stress = new tensor_t(e, 0);

    var.shpdx = new shapefn(e);
    if (NDIMS == 3) var.shpdy = new shapefn(e);
    var.shpdz = new shapefn(e);
}


static void create_matprops(const Param &par, Variables &var)
{
    var.mat = new MatProps(par, var);
}


void init(const Param& param, Variables& var)
{
    create_new_mesh(param, var);
    allocate_variables(var);
    create_matprops(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    *var.volume_old = *var.volume;
    compute_dvoldt(var, *var.tmp0, *var.dvoldt, *var.edvoldt);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    initial_stress_state(param, var, *var.stress, *var.strain, var.compensation_pressure);
    initial_weak_zone(param, var, *var.plstrain);
    initial_temperature(param, var, *var.temperature);
    apply_vbcs(param, var, *var.vel);
}


void update_mesh(const Param& param, Variables& var)
{
    update_coordinate(var, *var.coord);

    var.volume->swap(*var.volume_old);
    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    compute_dvoldt(var, *var.tmp0, *var.dvoldt, *var.edvoldt);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);
}


int main(int argc, const char* argv[])
{
    double start_time = 0;
#ifdef USE_OMP
    start_time = omp_get_wtime();
#endif

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

    var.max_vbc_val = find_max_vbc(param.bc);

    if (! param.sim.is_restarting) {
        init(param, var);
        output(param, var, start_time);
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
            output(param, var, start_time);
            var.frame ++;
        }

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    return 0;
}
