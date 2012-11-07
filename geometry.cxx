#include <cmath>
#include <limits>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "geometry.hpp"
#include "matprops.hpp"


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


void compute_volume(const double2d &coord, const int2d &connectivity,
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


double compute_dt(const Param& param, const Variables& var)
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

    double dt_elastic = 0.5 * minl / (param.bc.max_vbc_val * param.inertial_scaling);

    // std::cout << "dt: " << dt_maxwell << " " << dt_diffusion
    //           << " " << dt_elastic << "\n";

    return std::min(std::min(dt_elastic, dt_maxwell),
		    dt_diffusion);
}


void compute_mass(const Param &param,
                  const double2d &coord, const int2d &connectivity,
                  const double_vec &volume, const MatProps &mat,
                  double_vec &mass, double_vec &tmass)
{
    double pseudo_speed = param.bc.max_vbc_val * param.inertial_scaling;
    const int nelem = connectivity.shape()[0];
    for (int e=0; e<nelem; ++e) {
        double pseudo_rho = mat.bulkm(e) / (pseudo_speed * pseudo_speed);
        double m = pseudo_rho * volume[e] / NODES_PER_ELEM;
        double tm = mat.density(e) * mat.cp(e) * volume[e] / NODES_PER_ELEM;
        const int *conn = &connectivity[e][0];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            mass[conn[i]] += m;
            tmass[conn[i]] += tm;
        }
    }
    //for (int i=0; i<mass.size(); ++i)
    //    std::cout << i << ": mass = " << mass[i] << '\n';

    //for (int i=0; i<tmass.size(); ++i)
    //    std::cout << i << ": tmass = " << tmass[i] << '\n';
}


void compute_shape_fn(const double2d &coord, const int2d &connectivity,
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


