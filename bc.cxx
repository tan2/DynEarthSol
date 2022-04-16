#ifdef USE_NPROF
#include <nvToolsExt.h> 
#endif
#include <iostream>
#include <unordered_map>
#include <iomanip>
#include <math.h>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "markerset.hpp"
#include "barycentric-fn.hpp"
#include "geometry.hpp"
#include "utils.hpp"


#include "bc.hpp"


namespace {

//#pragma acc routine seq
void normal_vector_of_facet(int f, const int *conn, const array_t &coord,
                            double *normal, double &zcenter)
{
    int n0 = conn[NODE_OF_FACET[f][0]];
    int n1 = conn[NODE_OF_FACET[f][1]];
#ifdef THREED
    int n2 = conn[NODE_OF_FACET[f][2]];

    // vectors: n0-n1 and n0-n2
    double v01[NDIMS], v02[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = coord[n1][i] - coord[n0][i];
        v02[i] = coord[n2][i] - coord[n0][i];
    }

    // the outward normal vector is parallel to the cross product
    // of v01 and v02, with magnitude = area of the triangle
    normal[0] = (v01[1] * v02[2] - v01[2] * v02[1]) / 2;
    normal[1] = (v01[2] * v02[0] - v01[0] * v02[2]) / 2;
    normal[2] = (v01[0] * v02[1] - v01[1] * v02[0]) / 2;

    zcenter = (coord[n0][2] + coord[n1][2] + coord[n2][2]) / NODES_PER_FACET;
#else
    // the edge vector
    double v01[NDIMS];
    for (int i=0; i<NDIMS; ++i) {
        v01[i] = coord[n1][i] - coord[n0][i];
    }

    // the normal vector to the edge, pointing outward
    normal[0] = v01[1];
    normal[1] = -v01[0];

    zcenter = (coord[n0][1] + coord[n1][1]) / NODES_PER_FACET;
#endif
}

}

#pragma acc routine seq
bool is_on_boundary(const uint_vec &bcflag, int node)
{
    uint flag = bcflag[node];
    return flag & BOUND_ANY;
}


double find_max_vbc(const BC &bc, const double_vec &vbc_period_ratio_x)
{
    double max_vbc_val = 0;
    if (bc.vbc_x0 % 2 == 1) // odd number indicates fixed velocity component
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x0));
    if (bc.vbc_x1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x1));
    if (bc.vbc_y0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y0));
    if (bc.vbc_y1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y1));
    if (bc.vbc_z0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z0));
    if (bc.vbc_z1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z1));
    if (bc.vbc_n0 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n0));
    if (bc.vbc_n1 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n1));
    if (bc.vbc_n2 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n2));
    if (bc.vbc_n3 % 2 == 1)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_n3));

    return max_vbc_val;
}


void create_boundary_normals(const Variables &var, array_t &bnormals,
                             std::map<std::pair<int,int>, double*>  &edge_vectors)
{
    /* This subroutine finds the outward normal unit vectors of boundaries.
     * There are two types of boundaries: ibound{x,y,z}? and iboundn?.
     * The first type ibound{x,y,z}? has well defined "normal" direction
     * as their name implied. Thus they can be curved. (Ex. the top boundary
     * can have topography.)
     * The second type iboundn?, on the other hand, must be a plane so that
     * their vbc are well defined.
     *
     * If the normal component of the boundary velocity is fixed, the boundary
     * normal will not change with time.
     */

    for (int i=0; i<nbdrytypes; i++) {
        double normal[NDIMS] = {0};
        if (var.bfacets[i]->size() == 0) continue;

        for (auto j=var.bfacets[i]->begin(); j<var.bfacets[i]->end(); ++j) {
            int e = j->first;
            int f = j->second;
            double tmp;
            normal_vector_of_facet(f, (*var.connectivity)[e], *var.coord,
                                   normal, tmp);
            // make an unit vector
            double len = 0;
            for(int d=0; d<NDIMS; d++)
                len += normal[d] * normal[d];

            len = std::sqrt(len);
            for(int d=0; d<NDIMS; d++)
                normal[d] = normal[d] / len;

            if (j == var.bfacets[i]->begin()) {
                for(int d=0; d<NDIMS; d++)
                    bnormals[i][d] = normal[d];

                if (i < iboundn0) break; // slant boundaries start at iboundn0, other boundary can be curved
            }
            else {
                // Make sure the boundary is a plane, ie. all facets have the same normal vector.
                const double eps2 = 1e-12;
                double diff2 = 0;
                for(int d=0; d<NDIMS; d++)
                    diff2 += (bnormals[i][d] - normal[d]) * (bnormals[i][d] - normal[d]);
                if (diff2 > eps2) {
                    std::cerr << "Error: boundary " << i << " is curved.\n";
                    std::cerr << "Got  -  Expected\n";
                    std::cerr << bnormals[i][0] << " - " << normal[0] << '\n';
                    std::cerr << bnormals[i][1] << " - " << normal[1] << '\n';
#ifdef THREED
                    std::cerr << bnormals[i][2] << " - " << normal[2] << '\n';
#endif
                    std::exit(1);
                }
            }
        }
    }

    for (int i=0; i<nbdrytypes; i++) {
        if (var.bfacets[i]->size() == 0) continue;

        const double eps = 1e-15;
        for (int j=i+1; j<nbdrytypes; j++) {
            if (var.bfacets[j]->size() == 0) continue;
            double *s = new double[NDIMS];  // intersection of two boundaries
                                            // whole-application lifetime, no need to delete manually
#ifdef THREED
            // quick path: both walls are vertical
            if (std::abs((*var.bnormals)[i][NDIMS-1]) < eps &&
                std::abs((*var.bnormals)[j][NDIMS-1]) < eps) {
                s[0] = s[1] = 0;
                s[NDIMS-1] = 1;
            }
            else {
                // cross product of 2 normal vectors
                s[0] = (*var.bnormals)[i][1]*(*var.bnormals)[j][2] - (*var.bnormals)[i][2]*(*var.bnormals)[j][1];
                s[1] = (*var.bnormals)[i][2]*(*var.bnormals)[j][0] - (*var.bnormals)[i][0]*(*var.bnormals)[j][2];
                s[2] = (*var.bnormals)[i][0]*(*var.bnormals)[j][1] - (*var.bnormals)[i][1]*(*var.bnormals)[j][0];
            }
#else
            // 2D
            s[0] = 0;
            s[1] = 1;
#endif
            edge_vectors[std::make_pair(i, j)] = s;
        }
    }
}


void apply_vbcs(const Param &param, const Variables &var, array_t &vel)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // meaning of vbc flags (odd: free; even: fixed) --
    // 0: all components free
    // 1: normal component fixed, shear components free
    // 2: normal component free, shear components fixed at 0
    // 3: normal component fixed, shear components fixed at 0
    // 4: normal component free, shear component (not z) fixed, only in 3D
    // 5: normal component fixed at 0, shear component (not z) fixed, only in 3D

    const BC &bc = param.bc;
#ifdef THREED
    // vel period is not ready for 3D

#else

    // new way to find vbc_applied_x
    double t_now = var.time / YEAR2SEC;
    double vbc_applied_x0 = bc.vbc_val_x0 * interp1(bc.vbc_period_x0_time_in_yr,bc.vbc_period_x0_ratio, t_now);
    double vbc_applied_x1 = bc.vbc_val_x1 * interp1(bc.vbc_period_x1_time_in_yr,bc.vbc_period_x1_ratio, t_now);

    // find max and min coordinate for BOUNDX0
    double BOUNDX0_max = 0.;
    double BOUNDX0_min = 0.;
    double BOUNDX1_max = 0.;
    double BOUNDX1_min = 0.;
    bool if_init0 = false;
    bool if_init1 = false;
    for (int i=0; i<var.nnode; ++i) {
        if (! is_on_boundary(*var.bcflag, i)) continue;

        uint flag = (*var.bcflag)[i];
        if (flag & BOUNDX0) {
            const double *x = (*var.coord)[i];
            if (!if_init0) {
               BOUNDX0_max = x[1];
               BOUNDX0_min = x[1];
               if_init0 = true;
            } else {
                if (x[1]>BOUNDX0_max) {
                    BOUNDX0_max = x[1];
                }
                if (x[1]<BOUNDX0_min) {
                    BOUNDX0_min = x[1];
                }
            }
        } else if (flag & BOUNDX1) {
            const double *x = (*var.coord)[i];
            if (!if_init1) {
               BOUNDX1_max = x[1];
               BOUNDX1_min = x[1];
               if_init1 = true;
            } else {
                if (x[1]>BOUNDX1_max) {
                    BOUNDX1_max = x[1];
                }
                if (x[1]<BOUNDX1_min) {
                    BOUNDX1_min = x[1];
                }
            }            
        }
    }
    double BOUNDX0_width = BOUNDX0_max - BOUNDX0_min;
    double BOUNDX1_width = BOUNDX1_max - BOUNDX1_min;

    const double_vec& vbc_vertical_ratios_x0 = var.vbc_vertical_ratio_x0;
    const double_vec& vbc_vertical_ratios_x1 = var.vbc_vertical_ratio_x1;

    double_vec vbc_vertical_divisions_x0(4,0.);
    double_vec vbc_vertical_divisions_x1(4,0.);
    for (int i=0;i<4;i++) {
        vbc_vertical_divisions_x0[i] = - (BOUNDX0_max - var.vbc_vertical_div_x0[i] * BOUNDX0_width);
        vbc_vertical_divisions_x1[i] = - (BOUNDX0_max - var.vbc_vertical_div_x1[i] * BOUNDX0_width);
    }
#endif

    // diverging x-boundary
    const int var_nnode = var.nnode;
    const uint_vec *var_bcflag = var.bcflag;
    const array_t *var_coord = var.coord;
    const array_t *var_bnormals = var.bnormals;
    const int *var_vbc_types = var.vbc_types;
    const double *var_vbc_values = var.vbc_values;
    const std::map<std::pair<int,int>, double*>  *edgevec = &(var.edge_vectors);


    int bc_x0 = bc.vbc_x0;
    int bc_x1 = bc.vbc_x1;
    int bc_y0 = bc.vbc_y0;
    int bc_y1 = bc.vbc_y1;
    int bc_z0 = bc.vbc_z0;
    int bc_z1 = bc.vbc_z1;

    double bc_vx0 = bc.vbc_val_x0;
    double bc_vx1 = bc.vbc_val_x1;
    double bc_vy0 = bc.vbc_val_y0;
    double bc_vy1 = bc.vbc_val_y1;
    double bc_vz0 = bc.vbc_val_z0;
    double bc_vz1 = bc.vbc_val_z1;

    double bc_vx0min = bc.vbc_val_division_x0_min;
    double bc_vx0max = bc.vbc_val_division_x0_max;

    double bc_vx0r0 = bc.vbc_val_x0_ratio0;
    double bc_vx0r1 = bc.vbc_val_x0_ratio1;
    double bc_vx0r2 = bc.vbc_val_x0_ratio2;

#ifdef THREED
    #pragma omp parallel for default(none) \
        shared(bc, var, vel)
    #pragma acc parallel loop
#else
    #pragma omp parallel for default(none) \
        shared(bc, var, vel, vbc_vertical_ratios_x0, vbc_vertical_ratios_x1, \
                vbc_vertical_divisions_x0, vbc_vertical_divisions_x1, \
                vbc_applied_x0, vbc_applied_x1,var_bcflag,var_coord,var_bnormals,var_vbc_types, \
                var_vbc_values,edgevec,bc_x0,bc_vx0min,bc_vx0r0, \
                bc_vx0max,bc_vx0r2,bc_vx0r1,bc_x1,bc_z0,bc_z1,bc_vz0,bc_vz1)
    #pragma acc parallel loop
#endif
    for (int i=0; i<var_nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(*var_bcflag, i)) continue;

        uint flag = (*var_bcflag)[i];
        double *v = vel[i];

#ifdef THREED
#else
        const double *x = (*var_coord)[i];
        double ratio, rr, dvr;
#endif
        //
        // X
        //
        if (flag & BOUNDX0) {
            switch (bc_x0) {
            case 0:
                break;
            case 1:
#ifdef THREED
                v[0] = bc_vx0;
#else
                v[0] = vbc_applied_x0 * interp1(vbc_vertical_divisions_x0, vbc_vertical_ratios_x0, -x[1]);
#endif
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
#ifdef THREED
                v[0] = bc_vx0;
#else
                v[0] = vbc_applied_x0 * interp1(vbc_vertical_divisions_x0, vbc_vertical_ratios_x0,-x[1]);
#endif
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc_vx0;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc_vx0;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc_vx0;
                v[1] = 0;
                break;
#endif
            }
        }
        if (flag & BOUNDX1) {
            switch (bc_x1) {
            case 0:
                break;
            case 1:
#ifdef THREED
                v[0] = bc_vx1;
#else
                v[0] = vbc_applied_x1 * interp1(vbc_vertical_divisions_x1, vbc_vertical_ratios_x1, -x[1]);
#endif
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
#ifdef THREED
                v[0] = bc_vx1;
#else
                v[0] = vbc_applied_x1 * interp1(vbc_vertical_divisions_x1, vbc_vertical_ratios_x1, -x[1]);
#endif
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc_vx1;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc_vx1;
                v[2] = 0;
                break;
            case 7:
                v[0] = bc_vx1;
                v[1] = 0;
                break;
#endif
            }
        }
#ifdef THREED
        //
        // Y
        //
        if (flag & BOUNDY0) {
            switch (bc_y0) {
            case 0:
                break;
            case 1:
                v[1] = bc_vy0;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = 0;
                v[1] = bc_vy0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc_vy0;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc_vy0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc_vy0;
                break;
            }
        }
        if (flag & BOUNDY1) {
            switch (bc_y1) {
            case 0:
                break;
            case 1:
                v[1] = bc_vy1;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[0] = bc_vy1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc_vy1;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc_vy1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc_vy1;
                break;
            }
        }
#endif

        //
        // N
        //
        for (int ib=iboundn0; ib<=iboundn3; ib++) {
            const double eps = 1e-15;
            const double *n = (*var_bnormals)[ib]; // unit normal vector

            if (flag & (1 << ib)) {
                double fac = 0;
                switch (var_vbc_types[ib]) {
                case 1:
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS; d++)
                            vn += v[d] * n[d];  // normal velocity

			for (int d=0; d<NDIMS; d++)
                            v[d] += (var_vbc_values[ib] - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var_vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        vn += v[d] * n[d];  // normal velocity

				    for (int d=0; d<NDIMS; d++)
                                        v[d] += (var_vbc_values[ib] - vn) * n[d];  // setting normal velocity
                                }
                                else if (var_vbc_types[ic] == 1) {
                                    auto edge = edgevec->at(std::make_pair(ic, ib));
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 3:
                    for (int d=0; d<NDIMS; d++)
                        v[d] = var_vbc_values[ib] * n[d];  // v must be normal to n
                    break;
                case 11:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS-1; d++)
                            vn += v[d] * n[d];  // normal velocity

			for (int d=0; d<NDIMS-1; d++)
                            v[d] += (var_vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var_vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS-1; d++)
                                        vn += v[d] * n[d];  // normal velocity

				    for (int d=0; d<NDIMS-1; d++)
                                        v[d] += (var_vbc_values[ib] * fac - vn) * n[d];  // setting normal velocity
                                }
                                else if (var_vbc_types[ic] == 1) {
                                    auto edge = edgevec->at(std::make_pair(ic, ib));
                                    double ve = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        ve += v[d] * edge[d];

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] = ve * edge[d];  // v must be parallel to edge
                                }
                            }
                        }
                    }
                    break;
                case 13:
                    fac = 1 / std::sqrt(1 - n[NDIMS-1]*n[NDIMS-1]);  // factor for horizontal normal unit vector
                    for (int d=0; d<NDIMS-1; d++)
                        v[d] = var_vbc_values[ib] * fac * n[d];
                    v[NDIMS-1] = 0;
                    break;
                }
            }
        }

        //
        // Z, must be dealt last
        //

        // fast path: vz is usually free in the models
        if (bc_z0==0 && bc_z1==0) continue;

        if (flag & BOUNDZ0) {
            switch (bc_z0) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc_vz0;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc_vz0;
                break;
            }
        }
        if (flag & BOUNDZ1) {
            switch (bc_z1) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc_vz1;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc_vz1;
                break;
            }
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void apply_stress_bcs(const Param& param, const Variables& var, array_t& force)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // TODO: add general stress (Neumann) bcs

    if (param.control.gravity == 0) return;

    //
    // Gravity-induced (hydrostatic and lithostatic) stress BCs
    //
    for (int i=0; i<nbdrytypes; i++) {
        if (var.vbc_types[i] != 0 &&
            var.vbc_types[i] != 2 &&
            var.vbc_types[i] != 4) continue;

        if (i==iboundz0 && !param.bc.has_winkler_foundation) continue;
        if (i==iboundz1 && !param.bc.has_water_loading) continue;

        const auto& bdry = *(var.bfacets[i]);
        const auto& coord = *var.coord;

        const int bound = static_cast<int>(bdry.size());
        const conn_t *var_connectivity = var.connectivity;
        const array_t *var_coord = var.coord;
        const double var_compensation_pressure = var.compensation_pressure;
        const MatProps *var_mat = var.mat;
        const bool has_winkler_foundation = param.bc.has_winkler_foundation;
        const double winkler_delta_rho = param.bc.winkler_delta_rho;
        const double gravity = param.control.gravity;
        const double zlength = param.mesh.zlength;
        const double has_water_loading = param.bc.has_water_loading;
        const double surf_base_level = param.control.surf_base_level;

        // loops over all bdry facets
        #pragma acc parallel loop
        for (int n=0; n<bound; ++n) {
            // this facet belongs to element e
            int e = bdry[n].first;
            // this facet is the f-th facet of e
            int f = bdry[n].second;
            const int *conn = (*var_connectivity)[e];

            // the outward-normal vector
            double normal[NDIMS];
            // the z-coordinate of the facet center
            double zcenter;

            normal_vector_of_facet(f, conn, *var_coord, normal, zcenter);

            double p;
            if (i==iboundz0 && has_winkler_foundation) {
                // Winkler foundation for the bottom boundary
                p = var_compensation_pressure -
                    (var_mat->rho(e) + winkler_delta_rho) *
                    gravity * (zcenter + zlength);
            }
            else if (i==iboundz1 && has_water_loading) {
                // hydrostatic water loading for the surface boundary
                p = 0;
                if (zcenter < surf_base_level) {
                    // below sea level
                    const double sea_water_density = 1030;
                    p = sea_water_density * gravity * (surf_base_level - zcenter);
                }
            }
            else {
                // sidewalls
                p = ref_pressure(param, zcenter);
            }

            // lithostatc support - Archimed force (normal to the surface)
            for (int j=0; j<NODES_PER_FACET; ++j) {
                int nn = conn[NODE_OF_FACET[f][j]];
                double *f = force[nn];
                for (int d=0; d<NDIMS; ++d) {
                    #pragma acc atomic update
                    f[d] -= p * normal[d] / NODES_PER_FACET;
                }
            }
        }
    }

    if (param.bc.has_elastic_foundation) {
        /* A restoration force on the bottom nodes proportional to total vertical displacement */
        for (auto i=var.bnodes[iboundz0]->begin(); i<var.bnodes[iboundz0]->end(); ++i) {
            int n = *i;
            force[n][NDIMS-1] -= param.bc.elastic_foundation_constant * ((*var.coord)[n][NDIMS-1] - (*var.coord0)[n][NDIMS-1]);
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

    void simple_diffusion(const Variables& var)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        /* Diffusing surface topography to simulate the effect of erosion and
         * sedimentation.
         */

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;

        const int top_bdry = iboundz1;
        const auto& top = *(var.bfacets[top_bdry]);

        const int var_nnode = var.nnode;
        const int ntop = surfinfo.ntop;

        const conn_t *var_connectivity = var.connectivity;
#ifdef USE_NPROF
        nvtxRangePushA("prepare variable total_dx");
#endif
        double *total_dx = var.surfinfo.total_dx->data();
#ifdef USE_NPROF
        nvtxRangePop();
#endif
        double *total_slope = var.surfinfo.total_slope->data();
        double_vec& dh = *var.surfinfo.dh;

        #pragma acc parallel loop
        for (int i=0;i<var_nnode;i++) {
            total_dx[i] = 0.;
            total_slope[i] = 0.;
        }

        // loops over all top facets
#ifdef THREED
        const size_t tsize = top.size();
        #pragma acc parallel loop
        for (std::size_t i=0; i<tsize; ++i) {
            // this facet belongs to element e
            int e = top[i].first;
            // this facet is the f-th facet of e
            int f = top[i].second;

            const int *conn = (*var_connectivity)[e];
            int n0 = (*var_connectivity)[e][NODE_OF_FACET[f][0]];
            int n1 = (*var_connectivity)[e][NODE_OF_FACET[f][1]];

//#ifdef THREED
            int n2 = (*var_connectivity)[e][NODE_OF_FACET[f][2]];

            double projected_area;
            {
                // normal vector of this facet
                double normal[NDIMS];

                // two vectors n0-n1 and n0-n2
                // n is the cross product of these two vectors
                // the length of n is 2 * triangle area
                double x01, y01, z01, x02, y02, z02;
                x01 = coord[n1][0] - coord[n0][0];
                y01 = coord[n1][1] - coord[n0][1];
                z01 = coord[n1][2] - coord[n0][2];
                x02 = coord[n2][0] - coord[n0][0];
                y02 = coord[n2][1] - coord[n0][1];
                z02 = coord[n2][2] - coord[n0][2];

                normal[0] = y01*z02 - z01*y02;
                normal[1] = z01*x02 - x01*z02;
                normal[2] = x01*y02 - y01*x02;

                /* the area of this facet is:
                 *   0.5 * std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
                 *
                 * theta is the angle between this facet and horizontal
                 *   tan_theta = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]) / normal[2]
                 *   cos_theta = normal[2] / (2 * area)
                 *
                 * the area projected on the horizontal plane is:
                 *   projected_area = area * cos_theta = 0.5 * normal[2]
                 */
                projected_area = 0.5 * normal[2];
            }

            #pragma acc atomic update
            total_dx[n0] += projected_area;
            #pragma acc atomic update
            total_dx[n1] += projected_area;
            #pragma acc atomic update
            total_dx[n2] += projected_area;

            double shp2dx[NODES_PER_FACET], shp2dy[NODES_PER_FACET];
            double iv = 1 / (2 * projected_area);
            shp2dx[0] = iv * (coord[n1][1] - coord[n2][1]);
            shp2dx[1] = iv * (coord[n2][1] - coord[n0][1]);
            shp2dx[2] = iv * (coord[n0][1] - coord[n1][1]);
            shp2dy[0] = iv * (coord[n2][0] - coord[n1][0]);
            shp2dy[1] = iv * (coord[n0][0] - coord[n2][0]);
            shp2dy[2] = iv * (coord[n1][0] - coord[n0][0]);

            double D[NODES_PER_FACET][NODES_PER_FACET];
            for (int j=0; j<NODES_PER_FACET; j++) {
                for (int k=0; k<NODES_PER_FACET; k++) {
                    D[j][k] = (shp2dx[j] * shp2dx[k] +
                               shp2dy[j] * shp2dy[k]);
                }
            }

            const int n[NODES_PER_FACET] = {n0, n1, n2};
            for (int j=0; j<NODES_PER_FACET; j++) {
                double slope = 0;
                for (int k=0; k<NODES_PER_FACET; k++)
                    slope += D[j][k] * coord[n[k]][2];

                #pragma acc atomic update
                total_slope[n[j]] += slope * projected_area;
            }

            // std::cout << i << ' ' << n0 << ' ' << n1 << ' ' << n2 << "  "
            //           << projected_area << "  " << slope << '\n';
#else
            /* The 1D diffusion operation is implemented ad hoc, not using FEM
             * formulation (e.g. computing shape function derivation on the edges).
             */
        const size_t tsize = ntop-1;
        #pragma acc parallel loop
        for (std::size_t i=0; i<tsize;i++) {
            int n0 = top_nodes[i];
            int n1 = top_nodes[i+1];

            double dx = std::fabs(coord[n1][0] - coord[n0][0]);
            #pragma acc atomic update
            total_dx[n0] += dx;
            #pragma acc atomic update
            total_dx[n1] += dx;

            double slope = (coord[n1][1] - coord[n0][1]) / dx;
            #pragma acc atomic update
            total_slope[n0] -= slope;
            #pragma acc atomic update
            total_slope[n1] += slope;

            // std::cout << i << ' ' << n0 << ' ' << n1 << "  " << dx << "  " << slope << '\n';
#endif
        }

        const double var_dt = var.dt;
        const double surf_diff = surfinfo.surf_diff;
        const double base_level = surfinfo.base_level;
        const double diff_ratio_terrig = surfinfo.diff_ratio_terrig;
        const double diff_ratio_marine = surfinfo.diff_ratio_marine;
#ifdef THREED
        #pragma acc parallel loop
#endif
        for (int i=0; i<ntop; ++i) {
            // we don't treat edge nodes specially, i.e. reflecting bc is used for erosion.
            int n = top_nodes[i];
            double conv =  surf_diff * var_dt * total_slope[n] / total_dx[n];
#ifdef THREED
            dh[i] -= conv;
#else
            if ( coord[n][1] >  base_level && conv > 0.) {
                dh[i] -= diff_ratio_terrig * conv;
            } else if ( coord[n][1] < base_level && conv < 0. ) {
                dh[i] -= diff_ratio_marine * conv;
            } else {
                dh[i] -= conv;
            }
#endif
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void custom_surface_processes(const Variables& var, array_t& coord)
    {
        const int top_bdry = iboundz1;
        const int_vec& top_nodes = *var.bnodes[top_bdry];
        const std::size_t ntop = top_nodes.size();

        // loops over all top nodes
        for (std::size_t i=0; i<ntop; ++i) {
            int n = top_nodes[i];

            // coordinate of this node
            double x = coord[n][0];
            double z = coord[n][NDIMS-1];

            // compute topography changes here ...
            // positive dh means sedimentation, negative dh means erosion
            double dh = 0;
            {
                // dh == ....
            }

            coord[n][NDIMS-1] += dh;
        }
    }

    void get_surface_info(const Variables& var, \
        double_vec& top_base, double_vec& top_depth) {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();

        for (std::size_t i=0; i<ntop-1;i++) {
            int n0 = top_nodes[i];
            int n1 = top_nodes[i+1];

            double dx = std::fabs(coord[n1][0] - coord[n0][0]);
            top_base[i] += dx;
            top_base[i+1] += dx;
        }

        for (std::size_t i=0; i<ntop; i++) {
            top_depth[i] = surfinfo.base_level - coord[top_nodes[i]][1];
            top_base[i] *= 0.5;
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

    void out_basin_info(const int_vec& if_land) {
        const int chart_width = 80;
        const int data_num = if_land.size();
        const double ratio = double(chart_width) / data_num;

        int sum = 0;
        for (int i=0;i<data_num;i++) {
            double prec = i * ratio;
            if(prec > sum) {
                int diff = ceil(prec) - sum;
                sum += diff;
                for (int j=0;j<diff;j++)
                    std::cout << std::setw(1) << std::dec << if_land[i];
//                    printf("%d",if_land[i]);
            }
        }
        std::cout << " (ntop: " << std::setw(3) << data_num << ")\n";
//        printf(" (ntop: %d)\n",data_num);
    }

    void get_basin_info(const Variables& var, double_vec& top_depth, \
        std::vector<bool>& if_source, int_vec& if_land,\
        int_vec& starts, int_vec& ends, double_vec& dhacc_tmp, double_vec& dx) {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();
        int_vec left_mouse, right_mouse;
        int_vec starts0(2,0);
        left_mouse.reserve(4);
        right_mouse.reserve(4);

        starts0[0] = starts[0];
        starts0[1] = starts[1];

        // find land:1 and sea:0
        for (std::size_t i=0; i<ntop; i++) {
            double topo_tmp = top_depth[i] - dhacc_tmp[i];
            if (topo_tmp <= 0.) if_land[i] = 1;
        }

        // find river mouth
        for (std::size_t i=0; i<ntop-1; i++)
            if (if_land[i] + if_land[i+1] == 1) {
                if (if_land[i] == 1)
                    left_mouse.push_back(i);
                else
                    right_mouse.push_back(i+1);
            }
        
        // find basin sets
        int left_mouth_num = left_mouse.size();
        int right_mouth_num = right_mouse.size();

        if (left_mouth_num == 0) {
            if_source[0] = false;
        } else {
            if_source[0] = true;
            starts[0] = left_mouse[0];


            if (left_mouth_num == right_mouth_num) {
                if (left_mouse[0] < right_mouse[0]) {
                    ends[0] = right_mouse[0];
                } else {
                    if (left_mouth_num > 1) {
                        ends[0] = right_mouse[1];
                    } else {
                        ends[0] = ntop-3;
                    }
                }
            } else if (left_mouth_num < right_mouth_num) {
                ends[0] = right_mouse[1];
            } else {
                if (left_mouth_num == 1) {
                    ends[0] = ntop-3;
                } else {
                    ends[0] = right_mouse[0];
                }
            }
        }

        if (right_mouth_num == 0 ) {
            if_source[1] = false;
        } else {
            if_source[1] = true;
            starts[1] = right_mouse[right_mouth_num-1];

            if (left_mouth_num == right_mouth_num) {
                if (left_mouse[0] < right_mouse[0]) {
                    ends[1] = left_mouse[left_mouth_num-1];
                } else {
                    if (right_mouth_num > 1) {
                        ends[1] = left_mouse[left_mouth_num-2];
                    } else {
                        ends[1] = 2;
                    }
                }
            } else if (left_mouth_num > right_mouth_num) {
                ends[1] = left_mouse[left_mouth_num-2];
            } else {
                if (right_mouth_num == 1) {
                    ends[1] = 2;
                } else {
                    ends[1] = left_mouse[left_mouth_num-1];
                }
            }

        }
        for (int i=0;i<2;i++){
            if (if_source[i]) {
                int inext = pow(-1, i);
                double topo_tmp0 = top_depth[starts[i]] - dhacc_tmp[starts[i]];
                double topo_tmp1 = top_depth[starts[i]+inext] - dhacc_tmp[starts[i]+inext];
                double x0 = (*var.coord)[top_nodes[starts[i]]][0];
                double x1 = (*var.coord)[top_nodes[starts[i]+inext]][0];
                dx[i] = fabs( (x1-x0) * topo_tmp0 / (topo_tmp0 - topo_tmp1) );
//                printf("%d %d %f %f %f ",top_nodes[starts[i]],top_nodes[starts[i]+inext],x0,x1,dx[i]);
//                printf("%f %f \n",topo_tmp0,topo_tmp1);
            }
        }


        if ( var.steps%10000 == 0 ) {
            if (starts0[0] != starts[0] || starts0[1] != starts[1]) {
                printf("%d starts0: %d %d; starts: %d %d\n", \
                        var.steps,starts0[0],starts0[1],starts[0],starts[1]);
                out_basin_info(if_land);
//                for (std::size_t i=0;i<ntop;i++)
//                printf("%d",if_land[i]);
//                printf("\n");
            }
        }


        // the source is too close to boundary
        if (abs(starts[0]-int(ntop)/2) >= int(ntop)/2)
            if_source[0] = false;

        if (abs(starts[1]-int(ntop)/2) >= int(ntop)/2)
            if_source[1] = false;
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

    void simple_deposition(const Param& param,const Variables& var) {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
#ifdef THREED
        // not ready for 3D
        std::cout << "3D deposition of sediment processes is not ready yet.";
        exit(168);
//        std::cout << "Press enter to continue ...";
//        std::cin.get();
#endif
        double_vec &dh = *var.surfinfo.dh;
        double_vec &src_locs = *var.surfinfo.src_locs;
        double_vec &src_abj = *var.surfinfo.src_abj;

    #pragma acc serial
    for (int i=0;i<2;i++) {
        src_locs[i] = 0.;
        src_abj[i] = 0.;
    }

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();

        double_vec top_base(ntop,0.);
        double_vec top_depth(ntop,0.);

        std::vector<bool> if_source(2,true);
        double_vec dh_terrig(ntop,0.);
        double_vec dhacc_tmp(ntop,0.);
        int interval_report = 10000;

        if ( var.steps%interval_report == 0 )
            std::cout << "** Sedimentation report: ";
        get_surface_info(var,top_base,top_depth);

//******************************************************************
//              sedimentation by terrigenous source
//******************************************************************

        // The Pearl River Shibao et al. (2007)
        // 5.e7 ton/yr ~= 1. m^3/s
        // assuming the width is 50 km --> 2.e-5 m^2/s
        double terrig_width = 50.e3;
        double max_sedi_vol = param.control.surf_src_vol;
        if (max_sedi_vol != 1.)
            max_sedi_vol = max_sedi_vol / terrig_width * var.dt;
        else
            max_sedi_vol = param.control.surf_src_area * var.dt;

        int vol_ratio = 10;
        int ntry = 200;

        bool if_finish, if_shift;
        int nloop, nstep_in_basin;
        double unit_vol, depo_vol, dh_tmp, dist, old_dhacc_tmp;
        int_vec if_land(ntop,0);
        int_vec starts(2,0), ends(2,0);
        int_vec sign(2,0), nsedi(2,0);
        double_vec sedi_vol(2,0.), dx(2,0.);
        std::vector<bool> if_space_limited(2,false);
        std::vector<bool> if_slope_limited(2,false);

        sign[0] = 1;
        sign[1] = -1;

        get_basin_info(var,top_depth,if_source, if_land, starts, ends, dhacc_tmp, dx);

        // recored possible source locations
        for (int i=0;i<2;i++) {
//            printf("%d %d\n",if_source[0],if_source[1]);
            if (if_source[i]) {
                src_locs[i] = coord[top_nodes[starts[i]]][0] + dx[i] * pow(-1,i);
//                src_abj[i] = dx[i];
//                printf("side %d: %f\n",i, src_abj[i]);
            }
        }

        // deal with the sedimentation for the aspect of coastal source
        // to n times of both two sides
        for (int islide=0; islide<vol_ratio; islide++) {
            // do deposit for two sides
            for (int iside=0; iside<2; iside++) {
                // if iside is on source, do next side directly
                if (! if_source[iside]) continue;
                // of this side reach max_sedi_vol, do next side directly
                if (sedi_vol[iside] - max_sedi_vol > -1.e-4 ) continue;
                // slided volumn for each side
                unit_vol = std::min(max_sedi_vol / vol_ratio, max_sedi_vol - sedi_vol[iside]);
                // recorder of deposited volumn for do loop
                depo_vol = 0.;
                if_shift = true;
                // start do loop
                if_finish = false;
                nloop = 0;
                do {
                    // if coast line is changed or basin is filled, search topo for finding basin again
                    if (if_shift || if_space_limited[iside] || if_slope_limited[iside]) {
                        if_shift = false;
                        if_slope_limited[iside] = false;
                        dx[0] = 0.;
                        dx[1] = 0.;
                        get_basin_info(var,top_depth,if_source, if_land, starts, ends, dhacc_tmp, dx);
                        if (starts[iside] == ends[iside]) {
                            if_space_limited[iside] = true;
                            break;
                        } else {
                            if_space_limited[iside] = false;
                        }
                    }
                    // if this side become no source, do next side directly
                    if (! if_source[iside]) break;
                    // record times of one side loop
                    nloop++;
                    // the times for one side deposit
                    nstep_in_basin = abs(ends[iside] - starts[iside]);

                    // calculate the dh of basin
                    for (int istep=1; istep<nstep_in_basin; istep++) {

                        int j = starts[iside] + sign[iside]*istep;
                        int pj = starts[iside] + sign[iside]*(istep-1);
                        double slope_tmp = -1. * ( top_depth[j] - top_depth[j-sign[iside]] );
                        slope_tmp += dhacc_tmp[j] - dhacc_tmp[j-sign[iside]];

                        // stop propagate if slope is increasing
                        if ( slope_tmp > 0. ) {
                            if_slope_limited[iside] = true;
                            break;
                        }

                        // distance to the coastline
                        dist = fabs(coord[top_nodes[j]][0] - coord[top_nodes[starts[iside]]][0]); //Todo - dx[iside] + 1.;
                        if (dist < 1.) {
                            std::cout << dist << " ";
//                            printf("%f",dist);
                        }

//                        printf("%f\n",dist);
                        // deposit based on distance to coastline and remained unit volumn
                        dh_tmp = surfinfo.terrig_coeff * pow(surfinfo.terrig_base,nloop) * (unit_vol - depo_vol) * dist;

                        if ( dh_tmp != dh_tmp ) {
                            std::cout << "\ndh_tmp is NaN.\n";
                            std::cout << std::fixed << std::scientific << dh_tmp << "\t" << (unit_vol - depo_vol);
                            std::cout << "\t dist: " << std::dec << dist << std::endl;
                            std::cout << "\t j: " << j << std::endl;
                            std::cout << "\t top_nodes[j]: " << top_nodes[j] << std::endl;
                            std::cout << "\t coord[top_nodes[j]][0]: " << coord[top_nodes[j]][0] << std::endl;
                            std::cout << "\t iside: "<< iside << std::endl;
                            std::cout << "\t starts[iside]: " << starts[iside] << std::endl;
                            std::cout << "\ttop_nodes[starts[iside]]: " << top_nodes[starts[iside]] << std::endl;
                            std::cout << "\t coord[top_nodes[starts[iside]]][0]: " << coord[top_nodes[starts[iside]]][0] << std::endl;
                            exit(168);
//                            printf("\ndh_tmp is NaN.\n");
//                            printf("\n%e\t%e\t%f\n",dh_tmp, (unit_vol - depo_vol), dist);
                        }

                        old_dhacc_tmp = dhacc_tmp[j];
                        dhacc_tmp[j] += dh_tmp;

                        // if dhacc larger than depth, dhacc = depth (down is negative)
                        if (dhacc_tmp[j] > top_depth[j])
                            dhacc_tmp[j] = top_depth[j] + 0.1;


                        depo_vol += top_base[j] * (dhacc_tmp[j] - old_dhacc_tmp);

                        if (depo_vol >= unit_vol) {
                            if_finish = true;
                            break;
                        }
                    }

                    // move coast if basin is filled.
                    if (dhacc_tmp[starts[iside]+sign[iside]] >= top_depth[starts[iside]+sign[iside]]) {
                        starts[iside] += sign[iside];
                        if_shift = true;
                        if (starts[iside] == ends[iside])
                            if_space_limited[iside] = true;
                    }

                // end of do loop of one side deposit
                // finish of reach unit_vol or ntry
                } while (unit_vol - depo_vol > 1.e-4 && !if_finish && nloop <= ntry);

                nsedi[iside] += nloop;
                sedi_vol[iside] += depo_vol;
            } // end of each side

            if (if_space_limited[0] || if_space_limited[1]) break;

            for (std::size_t j=0; j<ntop; j++) {
                dh_terrig[j] += dhacc_tmp[j];
                dhacc_tmp[j] = 0.;
            }
            if (sedi_vol[0] >= max_sedi_vol && sedi_vol[1] >= max_sedi_vol) break;
        } // end of deposit slides

        if (if_source[0] || if_source[1]) {
            double_vec tmp_slope(ntop,0.);
            for (std::size_t i=0;i<ntop-1;i++) {
                tmp_slope[i] = 0.5 * ( dh_terrig[i+1] - dh_terrig[i] ) / top_base[i];
                tmp_slope[i+1] = 0.5 * ( dh_terrig[i+1] - dh_terrig[i] ) / top_base[i];
            }
            for (std::size_t i=0;i<ntop;i++)
                if (dh_terrig[i]!=0.)
                    if (fabs(tmp_slope[i]) > 3.e-3 ) {
                        std::cout << "arge ddh:" << std::setw(5) << std::dec << i << top_nodes[i];
                        std::cout <<  std::setw(9) << std::setprecision(2) << std::scientific << dh_terrig[i] << tmp_slope[i] << std::endl;
                    }
        }

        for (size_t i=0; i<ntop; i++)
            dh[i] += dh_terrig[i];

        if ( var.steps%interval_report == 0 ) {
            if (!if_source[0] || (!if_source[0]) ) {
                if (!if_source[0]) std::cout << "No source from 0. ";
                if (!if_source[1]) std::cout << "No source from 1. ";
                std::cout << std::endl;
            }
            for (int i=0;i<2;i++) {
                if (if_source[i]) {
                    if (if_space_limited[i]) std::cout << "\n   Space limited at " << i;
                    if (if_slope_limited[i]) std::cout << "\n   Slope limited at " << i;
                    std::cout << "\n    Side " << std::setw(1) << std::dec << i << ": Loc.: " << std::setw(8) << std::setprecision(2) << coord[top_nodes[starts[i]]][0]/1000.;
                    std::cout << " km (" << std::setw(5) << std::dec << starts[i] << " - " << std::setw(5) << ends[i] << "). Sediment: ";
                    std::cout << std::setw(8) << std::setprecision(2) << sedi_vol[i] << " m^2 (max: " << max_sedi_vol << " Loop: ";
                    std::cout << std::setw(5) << std::dec << nsedi[i] << std::endl;
                }
            }
            if (if_space_limited[0] || if_space_limited[1]) {
                std::cout << "\n    Space of basin is not enough for sediment . Do next round.";
                std::cout << std::setw(5) << std::dec << nsedi[0] << "/" << nsedi[1] << std::endl;
            }
        }

//******************************************************************
//              sedimentation by suspended source
//******************************************************************

        double ddh = surfinfo.depo_universal * var.dt;
        for (std::size_t i=0; i<ntop; i++)
            // if below the base level
            if (top_depth[i] > 0.) dh[i] += ddh;
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

    void simple_igneous(const Param& param,const Variables& var, double_vec& dh_oc, bool& has_partial_melting) {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
#ifdef THREED
        // not ready for 3D
        std::cout << "3D deposition of igneous processes is not ready yet.";
        exit(168);
#endif
        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();
        int_vec melting_marker;

        int nmarkers = var.markersets[0]->get_nmarkers();

        has_partial_melting = var.markersets[0]->if_melt(param.mat.mattype_partial_melting_mantle);

        if (!has_partial_melting) return;

        double_vec top_base(ntop,0.);
        double_vec top_depth(ntop,0.);

        get_surface_info(var,top_base,top_depth);

        double_vec melting_depth_elem(ntop,-101.e3);
        double_vec melting_depth(ntop, 0);

        for (int i=0;i<nmarkers;i++) {
            int m = var.markersets[0]->get_mattype(i);

            if (m == param.mat.mattype_partial_melting_mantle) {
                melting_marker.push_back(i); 

                int e = var.markersets[0]->get_elem(i);
                double vol = (*var.volume)[e];
                int num_markers_in_elem = 0;

                for( int k = 0; k < param.mat.nmat; k++ )
                    num_markers_in_elem += (*var.elemmarkers)[e][k];

                double x[NDIMS] = {0};
                for (int j = 0; j < NDIMS; j++) {
                    for (int k = 0; k < NODES_PER_ELEM; k++)
                        x[j] += var.markersets[0]->get_eta(i)[k]*
                            (*var.coord)[ (*var.connectivity)[e][k] ][j];
                }

                // search top element node for melting marker
                for (size_t j=1; j<ntop;j++) {
                    double dx0 = x[0] - coord[top_nodes[j-1]][0];
                    double dx1 = x[0] - coord[top_nodes[j]][0];
                    if (dx0*dx1 >= 0) continue;

                    double ratio = dx0/(dx0-dx1);
                    double dvol = param.mat.convert_rate_oceanic_crust * vol / num_markers_in_elem * var.dt;// * param.mesh.quality_check_step_interval;

                    dh_oc[j-1] += ratio * dvol / top_base[j-1];
                    dh_oc[j] += (1.-ratio) * dvol / top_base[j];
                    // find the melting tops on each element facet
                    if (melting_depth_elem[j-1] < x[NDIMS-1] ) {
                        melting_depth_elem[j-1] = x[NDIMS-1];
                    }
                }
            }
        }
        if (melting_marker.size() == 0)
            return;
        // find the depth of melting top of node
        for (size_t i=1;i<ntop-1;i++) {
            if (melting_depth_elem[i-1] < -100.e3 && melting_depth_elem[i] < -100.e3) continue;
            if (melting_depth_elem[i-1] < -100.e3)
                melting_depth[i] = melting_depth_elem[i];
            else if (melting_depth_elem[i] < -100.e3)
                melting_depth[i] = melting_depth_elem[i-1];
            else
                melting_depth[i] = 0.5 * ( melting_depth_elem[i-1] + melting_depth_elem[i] );
        }
        // find the max depth of melting top
        double max_dp = 0.;
        for (size_t i=0;i<ntop;i++)
            max_dp = std::min(max_dp, melting_depth[i]);
        // normalizing the depth of melting
        for (size_t i=0;i<ntop;i++) {
            if (melting_depth[i] == 0.)
                melting_depth[i] = max_dp;
            melting_depth_elem[i] = melting_depth[i] / max_dp;
        }
        // smoothing the depth of melting
        for (size_t i=1;i<ntop-1;i++)
            melting_depth[i] = (melting_depth_elem[i-1] + melting_depth_elem[i] + melting_depth_elem[i+1])/3.;
        melting_depth[0] = melting_depth[1];
        melting_depth[ntop-1] = melting_depth[ntop-2];
        for (size_t i=0;i<ntop;i++)
            dh_oc[i] = dh_oc[i] / melting_depth[i];
        if (var.steps % param.mesh.quality_check_step_interval == 0)
            if (melting_marker.size() > 0) {
                for (size_t i=0;i<ntop;i++)
                    if (dh_oc[i] > 0) {
                        double coef = dh_oc[i]/ var.dt * 1000. * YEAR2SEC; // / param.mesh.quality_check_step_interval ;
                        printf("%zu x: %f dh_oc: %f (mm/yr) depth: %f\n",i, coord[top_nodes[i]][0],coef, melting_depth[i]);
                    }
            }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }
}


void surface_plstrain_diffusion(const Param &param, \
    const Variables& var, double_vec& plstrain)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double half_life = 1.e2 * YEAR2SEC;
    double lambha = 0.69314718056 / half_life; // ln2
    #pragma omp parallel for default(none)      \
        shared(param, var, plstrain, lambha)
    for (auto e=(*var.top_elems).begin();e<(*var.top_elems).end();e++) {
        // Find the most abundant marker mattype in this element
        int_vec &a = (*var.elemmarkers)[*e];
        int mat = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
        if (mat != param.mat.mattype_oceanic_crust)
            plstrain[*e] -= plstrain[*e] * lambha * var.dt;
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void correct_surface_element(const Variables& var, \
    const double_vec& dhacc, MarkerSet& ms, tensor_t& stress, \
    tensor_t& strain, tensor_t& strain_rate, double_vec& plstrain)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    const size_t ntop_elem = var.top_elems->size();
    array_t coord0s(ntop_elem*NODES_PER_ELEM,0.);
    double_vec new_volumes(ntop_elem,0.);

#ifdef LLVM
    #pragma omp parallel for default(none)      \
        shared(ntop_elem, var, dhacc, stress,   \
        strain, strain_rate, plstrain, coord0s, new_volumes)
#else
    #pragma omp parallel for default(none)      \
        shared(var, dhacc, stress, strain,      \
        strain_rate, plstrain, coord0s, new_volumes)
#endif
    for (size_t i=0;i<ntop_elem;i++) {
        const double *coord1[NODES_PER_ELEM];

        auto e = (*var.top_elems)[i];
        int* tnodes = (*var.connectivity)[e];
        double *c00 = coord0s[i*NODES_PER_ELEM];
        double *c01 = coord0s[i*NODES_PER_ELEM+1];
        double *c02 = coord0s[i*NODES_PER_ELEM+2];

        for (int j=0; j<NODES_PER_ELEM;j++)
            coord1[j] = (*var.coord)[tnodes[j]];
        compute_volume(coord1, new_volumes[i]);

        // restore the reference node locations before deposition/erosion 
        c00[0] = (*var.coord)[tnodes[0]][0];
        c00[1] = (*var.coord)[tnodes[0]][1] - dhacc[tnodes[0]];
        c01[0] = (*var.coord)[tnodes[1]][0];
        c01[1] = (*var.coord)[tnodes[1]][1] - dhacc[tnodes[1]];
        c02[0] = (*var.coord)[tnodes[2]][0];
        c02[1] = (*var.coord)[tnodes[2]][1] - dhacc[tnodes[2]];

        // correct stress and strain
        double dArea0 = ((c01[0] - c00[0])*(c02[1] - c00[1]) \
                       - (c02[0] - c00[0])*(c01[1] - c00[1]))/2.0;
        dArea0 = (dArea0 > 0.0) ? dArea0 : -dArea0;

        double rdv = new_volumes[i] / dArea0;

        if (rdv > 1.) {
            plstrain[e] /= rdv;
            for (int i=0;i<NSTR;i++) {
                stress[e][i] /= rdv;
                strain[e][i] /= rdv;
                strain_rate[e][i] /= rdv;
            }
        }
    }
    Barycentric_transformation bary(*var.top_elems, *var.coord, *var.connectivity, new_volumes);
    ms.correct_surface_marker(var, coord0s, bary);
#ifdef USE_NPROF
    nvtxRangePop();
#endif

}

void surface_processes(const Param& param, const Variables& var, array_t& coord, tensor_t& stress, tensor_t& strain, \
                       tensor_t& strain_rate, double_vec& plstrain, SurfaceInfo& surfinfo, \
                        std::vector<MarkerSet*> &markersets, int_vec2D& elemmarkers)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    const int slow_updates_interval = 10;
    bool has_partial_melting = false;
    const int ntop = var.surfinfo.ntop;
    double_vec &dh = *var.surfinfo.dh;
    double_vec &dh_oc = *var.surfinfo.dh_oc;

    #pragma acc parallel loop
    for (int i=0;i<ntop;i++)
        dh[i] = 0.;

    switch (param.control.surface_process_option) {
    case 0:
        // no surface process
        break;
    case 1:
        simple_diffusion(var);
        break;
    case 101:
        custom_surface_processes(var, coord);
        break;
    case 102:
#ifdef THREED
        std::cout << "3D deposition of sediment processes is not ready yet.";
        exit(168);
#else
        simple_diffusion(var);
        if (var.steps != 0)
            simple_deposition(param, var);
#endif
        break;
    default:
        std::cout << "Error: unknown surface process option: " << param.control.surface_process_option << '\n';
        std::exit(1);
    }
#ifdef THREED
    const int_vec *top_nodes = surfinfo.top_nodes;
    #pragma acc parallel loop
    for (std::size_t i=0; i<ntop; ++i) {
        int n = (*top_nodes)[i];
        coord[n][NDIMS-1] += dh[i];
    }
    // todo
#else
    // go through all surface nodes and abject all surface node by dh
    for (int i=0; i<ntop; i++) {
        // get global index of node
        // update coordinate via dh
        // update dhacc for marker correction
        int n = (*surfinfo.top_nodes)[i];
        (coord)[n][NDIMS-1] += dh[i];
        (*surfinfo.dhacc)[n] += dh[i];

        // go through connected elements
        // get local index of surface element
        // get global index of element
        // get local index of node in connected element
        // update edhacc of connected elements
        for (std::size_t j=0; j<(*surfinfo.node_and_elems)[i].size(); j++) {
            int e = (*surfinfo.node_and_elems)[i][j];
            int eg = (*surfinfo.top_facet_elems)[e];
            int ind = (*surfinfo.arcelem_and_nodes_num)[e][i];
            (*surfinfo.edhacc)[eg][ind] += dh[i];
        }
    }
#endif

#ifdef THREED

#else
    if (var.steps != 0) {
        if ( var.steps % param.mesh.quality_check_step_interval == 0) {
            // correct surface marker.
            correct_surface_element(var, *surfinfo.dhacc, *markersets[0], stress, strain, strain_rate, plstrain);
            std::fill(surfinfo.dhacc->begin(), surfinfo.dhacc->end(), 0.);
            // set marker of sediment.
            markersets[0]->set_surface_marker(var, param.mat.mattype_sed, *surfinfo.edhacc, elemmarkers);
        }
    }
#endif
    if ( param.mat.phase_change_option == 2) {
#ifdef THREED
        std::cout << "3D simple_igneous processes is not ready yet.";
        exit(168);
#else
        simple_igneous(param,var, dh_oc, has_partial_melting);

        if ( has_partial_melting ) {
            // go through all surface nodes and abject all surface node by dh_oc
            // as same as dh
            for (int i=0; i<ntop; i++) {
                int n = (*surfinfo.top_nodes)[i];
                (coord)[n][NDIMS-1] += dh_oc[i];
                (*surfinfo.dhacc_oc)[n] += dh_oc[i];

                for (std::size_t j=0; j<(*surfinfo.node_and_elems)[i].size(); j++) {
                    int e = (*surfinfo.node_and_elems)[i][j];
                    int eg = (*surfinfo.top_facet_elems)[e];
                    int ind = (*surfinfo.arcelem_and_nodes_num)[e][i];
                    (*surfinfo.edhacc_oc)[eg][ind] += dh_oc[i];
                }
            }

            if (!(var.steps % param.mesh.quality_check_step_interval)) {
                // correct surface marker.
                correct_surface_element(var, *surfinfo.dhacc_oc, *markersets[0], stress, strain, strain_rate, plstrain);
                std::fill(surfinfo.dhacc_oc->begin(), surfinfo.dhacc_oc->end(), 0.);
                // set marker of sediment.
                markersets[0]->set_surface_marker(var,param.mat.mattype_oceanic_crust,*surfinfo.edhacc_oc,elemmarkers);
            }
        }
#endif
    }

#ifdef THREED
    // todo
#else
    if (!(var.steps % param.mesh.quality_check_step_interval &&  var.steps != 0))
        // correct the plastic strain of urface element for preventing surface landslide.
        surface_plstrain_diffusion(param,var,plstrain);

    if ( var.steps%10000 == 0 &&  var.steps != 0  ) {
        double max_dh = 0., min_dh = 0., max_dh_oc = 0.;
        for (int i=0;i<ntop;i++) {
            max_dh = std::max(max_dh, dh[i]);
            min_dh = std::min(min_dh, dh[i]);
            // std::cout << n << "  dh:  " << dh << '\n';
        }

        std::cout << "   max erosion / sedimentation rate (mm/yr):  "
                    << std::fixed << std::setprecision(3) << min_dh / var.dt * 1000. * YEAR2SEC << " / "
                    << std::fixed << std::setprecision(3) << max_dh / var.dt * 1000. * YEAR2SEC << '\n';

        if ( param.mat.phase_change_option == 2) {
            for (int i=0;i<ntop;i++)
                max_dh_oc = std::max(max_dh_oc, dh_oc[i]);
            std::cout << "   max igneous eruption rate (mm/yr):  "
                        << std::fixed << std::setprecision(3) << max_dh_oc / var.dt * 1000. * YEAR2SEC << '\n';
        }
    }
#endif

#ifdef USE_NPROF
    nvtxRangePop();
#endif


}
