#include <cstdio>
#include <ctime>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
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


static double tetrahedron_volume(const double *d0,
                                 const double *d1,
                                 const double *d2,
                                 const double *d3)
{
    double x01 = d0[0] - d1[0];
    double x12 = d1[0] - d2[0];
    double x23 = d2[0] - d3[0];

    double y01 = d0[1] - d1[1];
    double y12 = d1[1] - d2[1];
    double y23 = d2[1] - d3[1];

    double z01 = d0[2] - d1[2];
    double z12 = d1[2] - d2[2];
    double z23 = d2[2] - d3[2];

    return (x01*(y23*z12 - y12*z23) +
            x12*(y01*z23 - y23*z01) +
            x23*(y12*z01 - y01*z12)) / 6;
}


static double triangle_area(const double *a,
                            const double *b,
                            const double *c)
{
    double ab0, ab1, ac0, ac1;

    // ab: vector from a to b
    ab0 = b[0] - a[0];
    ab1 = b[1] - a[1];
    // ac: vector from a to c
    ac0 = c[0] - a[0];
    ac1 = c[1] - a[1];

    // area = (cross product of ab and ac) / 2
    return (ab0*ac1 - ab1*ac0) / 2;
}


static void compute_volume(const double2d &coord, const int2d &connectivity,
                           double_vec &volume, double_vec &volume_n)
{
    const int nelem = connectivity.shape()[0];
    for (int e=0; e<nelem; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = &coord[n0][0];
        const double *b = &coord[n1][0];
        const double *c = &coord[n2][0];

        double vol;
        if (NDIMS == 3) {
            int n3 = connectivity[e][3];
            const double *d = &coord[n3][0];
            vol = tetrahedron_volume(a, b, c, d);
        }
        else {
            vol = triangle_area(a, b, c);
        }
        volume[e] = vol;
        //std::cout << e << ": volume =" << vol << '\n';
    }

    // volume_n is (node-averaged volume * NODES_PER_ELEM)
    // volume_n[n] is init'd to 0 by resize()
    for (int e=0; e<nelem; ++e) {
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = connectivity[e][i];
            volume_n[n] += volume[e];
        }
    }

    //for (int i=0; i<volume_n.size(); ++i)
    //    std::cout << i << ": volume_n = " << volume_n[i] << '\n';
}


static void compute_mass(const double2d &coord, const int2d &connectivity,
                         const double_vec &volume, const MatProps &mat,
                         double_vec &mass, double_vec &tmass)
{
    const int nelem = connectivity.shape()[0];
    for (int e=0; e<nelem; ++e) {
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = connectivity[e][i];
            // TODO
            const double maxvbcval = 1e-10;
            const double pseudo_factor = 1e5; // == 1/strain_inert in geoflac
            double pseudo_speed = maxvbcval * pseudo_factor;
            double pseudo_rho = mat.bulkm(e) / (pseudo_speed * pseudo_speed);
            mass[n] += pseudo_rho * volume[e] / NODES_PER_ELEM;
            tmass[n] += mat.density(e) * mat.cp(e) * volume[e] / NODES_PER_ELEM;
        }
    }
    //for (int i=0; i<mass.size(); ++i)
    //    std::cout << i << ": mass = " << mass[i] << '\n';

    //for (int i=0; i<tmass.size(); ++i)
    //    std::cout << i << ": tmass = " << tmass[i] << '\n';
}


static void compute_shape_fn(const double2d &coord, const int2d &connectivity,
                             const double_vec &volume,
                             double2d &shpdx, double2d &shpdy, double2d &shpdz)
{
    const int nelem = connectivity.shape()[0];
    for (int e=0; e<nelem; ++e) {

        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *d0 = &coord[n0][0];
        const double *d1 = &coord[n1][0];
        const double *d2 = &coord[n2][0];

        if (NDIMS == 3) {
            int n3 = connectivity[e][3];
            const double *d3 = &coord[n3][0];

            double iv = 1 / (6 * volume[e]);

            double x01 = d0[0] - d1[0];
            double x02 = d0[0] - d2[0];
            double x03 = d0[0] - d3[0];
            double x12 = d1[0] - d2[0];
            double x13 = d1[0] - d3[0];
            double x23 = d2[0] - d3[0];

            double y01 = d0[1] - d1[1];
            double y02 = d0[1] - d2[1];
            double y03 = d0[1] - d3[1];
            double y12 = d1[1] - d2[1];
            double y13 = d1[1] - d3[1];
            double y23 = d2[1] - d3[1];

            double z01 = d0[2] - d1[2];
            double z02 = d0[2] - d2[2];
            double z03 = d0[2] - d3[2];
            double z12 = d1[2] - d2[2];
            double z13 = d1[2] - d3[2];
            double z23 = d2[2] - d3[2];

            shpdx[e][0] = iv * (y13*z12 - y12*z13);
            shpdx[e][1] = iv * (y02*z23 - y23*z02);
            shpdx[e][2] = iv * (y13*z03 - y03*z13);
            shpdx[e][3] = iv * (y01*z02 - y02*z01);

            shpdy[e][0] = iv * (z13*x12 - z12*x13);
            shpdy[e][1] = iv * (z02*x23 - z23*x02);
            shpdy[e][2] = iv * (z13*x03 - z03*x13);
            shpdy[e][3] = iv * (z01*x02 - z02*x01);

            shpdz[e][0] = iv * (x13*y12 - x12*y13);
            shpdz[e][1] = iv * (x02*y23 - x23*y02);
            shpdz[e][2] = iv * (x13*y03 - x03*y13);
            shpdz[e][3] = iv * (x01*y02 - x02*y01);
        }
        else {
            double iv = 1 / (2 * volume[e]);

            shpdx[e][0] = iv * (d1[1] - d2[1]);
            shpdx[e][1] = iv * (d2[1] - d0[1]);
            shpdx[e][2] = iv * (d0[1] - d1[1]);

            shpdz[e][0] = iv * (d2[0] - d1[0]);
            shpdz[e][1] = iv * (d0[0] - d2[0]);
            shpdz[e][2] = iv * (d1[0] - d0[0]);
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


void init(const Param& param, Variables& var)
{
    void create_matprops(const Param&, Variables&);

    create_new_mesh(param, var);
    allocate_variables(var);
    // XXX
    //create_nsupport(*connectivity, nnode, var.nsupport, var.support);
    create_matprops(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.volume, *var.volume_n);
    compute_mass(*var.coord, *var.connectivity, *var.volume, *var.mat,
                 *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume,
                     *var.shpdx, *var.shpdy, *var.shpdz);
    // XXX
    //create_jacobian();

    initial_temperature(param, var, *var.temperature);
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
void update_mesh() {};
void rotate_stress() {};


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
    f = fopen(buffer, "w");
    fwrite(var.coord->data(), sizeof(double), var.coord->num_elements(), f);
    fclose(f);

    // connectivity
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "connectivity", var.frame);
    f = fopen(buffer, "w");
    fwrite(var.connectivity->data(), sizeof(int), var.connectivity->num_elements(), f);
    fclose(f);

    // temperature
    snprintf(buffer, 255, "%s.%s.%06d", param.sim.modelname.c_str(), "temperature", var.frame);
    f = fopen(buffer, "w");
    fwrite(&(var.temperature->front()), sizeof(double), var.temperature->size(), f);
    fclose(f);

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
    var.dt = 1e7;
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

    do {
        var.steps ++;
        var.time += var.dt;

        update_temperature(param, var, *var.temperature, *var.tmp0);
        update_strain_rate();
        update_stress();
        update_force();
        update_mesh();
        rotate_stress();

        if ( (var.steps == var.frame * param.sim.output_step_interval) ||
             (var.time > var.frame * param.sim.output_time_interval_in_yr * YEAR2SEC) ) {
            output(param, var);
            std::cout << "  Output # " << var.frame
                      << " , step = " << var.steps
                      << " , time = " << var.time / YEAR2SEC << " yr.\n";
            var.frame ++;
        }

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    return 0;
}
