#include <cstdio>
#include <ctime>
#include <limits>
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

#ifndef THREED
    // area = norm(cross product of ab and ac) / 2
    return std::fabs(ab0*ac1 - ab1*ac0) / 2;
#else
    double ab2, ac2;
    ab2 = b[2] - a[2];
    ac2 = c[2] - a[2];

    // vector components of ab x ac
    double d0, d1, d2;
    d0 = ab1*ac2 - ab2*ac1;
    d1 = ab2*ac0 - ab0*ac2;
    d2 = ab0*ac1 - ab1*ac0;
    
    // area = norm(cross product of ab and ac) / 2
    return std::sqrt(d0*d0 + d1*d1 + d2*d2) / 2;
#endif
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


/* Give two points, return square of the distance */
static double dist2(const double* a, const double* b)
{
    double sum = 0;;
    for (int i=0; i<NDIMS; ++i) {
        double d = b[i] - a[i];
        sum += d * d;
    }
    return sum;
}


static double compute_dt(const Param& param, const Variables& var)
{
    const int nelem = var.nelem;
    const int2d_ref& connectivity = *var.connectivity;
    const double2d_ref& coord = *var.coord;
    const double_vec& volume = *var.volume;

    double dt_maxwell = std::numeric_limits<double>::max();
    double dt_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();

    for (int e=0; e<nelem; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = &coord[n0][0];
        const double *b = &coord[n1][0];
        const double *c = &coord[n2][0];

        // min height of this element
        double minh;
        if (NDIMS == 3) {
            int n3 = connectivity[e][3];
            const double *d = &coord[n3][0];

            // max facet area of this tet
            double maxa = std::max(std::max(triangle_area(a, b, c),
                                            triangle_area(a, b, d)),
                                   std::max(triangle_area(c, d, a),
                                            triangle_area(c, d, b)));
            minh = 3 * volume[e] / maxa;
        }
        else {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * volume[e] / maxl;

        }
        dt_maxwell = std::min(dt_maxwell,
                              0.5 * var.mat->visc_min / (1e-40 + var.mat->shearm(e)));
        dt_diffusion = std::min(dt_diffusion,
                                0.5 * minh * minh / var.mat->therm_diff_max);
	minl = std::min(minl, minh);
    }

    double dt_elastic = 0.5 * param.strain_inert * minl / param.maxvbcval;

    // std::cout << "dt: " << dt_maxwell << " " << dt_diffusion
    //           << " " << dt_elastic << "\n";

    return std::min(std::min(dt_elastic, dt_maxwell),
		    dt_diffusion);
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


void apply_vbcs(const Param &param, const Variables &var, double2d &vel)
{
    const double maxvbcval = 1e-10;
    // TODO: adding different types of vbcs later

    // diverging x-boundary
    for (int i=0; i<var.nnode; ++i) {
        int flag = (*var.bcflag)[i];

        // X
        if (flag & BOUNDX0) {
            vel[i][0] = -maxvbcval;
        }
        else if (flag & BOUNDX1) {
            vel[i][0] = maxvbcval;
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

    // XXX
    param.strain_inert = 1e-5;
    param.maxvbcval = 1e-10;
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
