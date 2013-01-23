#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "bc.hpp"


bool is_on_boundary(const Variables &var, int node)
{
    int flag = (*var.bcflag)[node];
    return flag & (BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | BOUNDZ0 | BOUNDZ1);
}


double find_max_vbc(const BC &bc)
{
    double max_vbc_val = 0;
    if (bc.vbc_x0 == 1 || bc.vbc_x0 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x0));
    if (bc.vbc_x1 == 1 || bc.vbc_x1 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x1));
    if (bc.vbc_y0 == 1 || bc.vbc_y0 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y0));
    if (bc.vbc_y1 == 1 || bc.vbc_y1 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y1));
    if (bc.vbc_z0 == 1 || bc.vbc_z0 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z0));
    if (bc.vbc_z1 == 1 || bc.vbc_z1 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z1));

    return max_vbc_val;
}


void apply_vbcs(const Param &param, const Variables &var, array_t &vel)
{
    // meaning of vbc flags (odd: free; even: fixed) --
    // 0: all components free
    // 1: normal component fixed, shear components free
    // 2: normal component free, shear components fixed at 0
    // 3: normal component fixed, shear components fixed at 0
    // 4: normal component free, shear component (not z) fixed, only in 3D
    // 5: normal component fixed at 0, shear component (not z) fixed, only in 3D

    const BC &bc = param.bc;

    // diverging x-boundary
    #pragma omp parallel for default(none) \
        shared(bc, var, vel)
    for (int i=0; i<var.nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(var, i)) continue;

        int flag = (*var.bcflag)[i];
        double *v = vel[i];

        // X
        if (flag & BOUNDX0) {
            switch (bc.vbc_x0) {
            case 0:
                break;
            case 1:
                v[0] = bc.vbc_val_x0;
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
                v[0] = bc.vbc_val_x0;
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc.vbc_val_x0;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc.vbc_val_x0;
                v[2] = 0;
                break;
#endif
            }
        }
        else if (flag & BOUNDX1) {
            switch (bc.vbc_x1) {
            case 0:
                break;
            case 1:
                v[0] = bc.vbc_val_x1;
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
                v[0] = bc.vbc_val_x1;
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc.vbc_val_x1;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc.vbc_val_x1;
                v[2] = 0;
                break;
#endif
            }
        }
#ifdef THREED
        // Y
        if (flag & BOUNDY0) {
            switch (bc.vbc_y0) {
            case 0:
                break;
            case 1:
                v[1] = bc.vbc_val_y0;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[1] = bc.vbc_val_y0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            case 5:
                v[1] = 0;
                v[0] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            }
        }
        else if (flag & BOUNDY1) {
            switch (bc.vbc_y1) {
            case 0:
                break;
            case 1:
                v[1] = bc.vbc_val_y1;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[1] = bc.vbc_val_y1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y1;
                v[2] = 0;
                break;
            case 5:
                v[1] = 0;
                v[0] = bc.vbc_val_y1;
                v[2] = 0;
                break;
            }
        }
#endif
        // Z

        // fast path: vz is usually free in the models
        if (bc.vbc_z0==0 && bc.vbc_z1==0) continue;

        if (flag & BOUNDZ0) {
            switch (bc.vbc_z0) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc.vbc_val_z0;
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
                v[NDIMS-1] = bc.vbc_val_z0;
                break;
            case 4:
                v[0] = bc.vbc_val_z0;
#ifdef THREED
                v[1] = bc.vbc_val_z0;
#endif
                break;
            case 5:
                v[0] = bc.vbc_val_z0;
#ifdef THREED
                v[1] = bc.vbc_val_z0;
#endif
                v[NDIMS-1] = 0;
                break;
            }
        }
        else if (flag & BOUNDZ1) {
            switch (bc.vbc_z1) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc.vbc_val_z1;
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
                v[NDIMS-1] = bc.vbc_val_z1;
                break;
            case 4:
                v[0] = bc.vbc_val_z1;
#ifdef THREED
                v[1] = bc.vbc_val_z1;
#endif
                break;
            case 5:
                v[0] = bc.vbc_val_z1;
#ifdef THREED
                v[1] = bc.vbc_val_z1;
#endif
                v[NDIMS-1] = 0;
                break;
            }
        }
    }
}


void apply_stress_bcs(const Param& param, const Variables& var, array_t& force)
{
    // TODO: add general stress (Neumann) bcs

    // hydrostatic water loading for the surface boundary
    if (param.bc.water_loading && param.control.gravity != 0) {
        std::cout << "TODO: applying water loading...\n";
    }

    // Wrinkler foundation for the bottom boundary
    if (param.bc.wrinkler_foundation && param.control.gravity != 0) {
        const int bottom_bdry = bdry_order.find(BOUNDZ0)->second;
        const auto& bottom = var.bfacets[bottom_bdry];
        const auto& coord = *var.coord;
        // loops over all bottom facets
        for (std::size_t i=0; i<bottom.size(); ++i) {
            // this facet belongs to element e
            int e = bottom[i].first;
            // this facet is the f-th facet of e
            int f = bottom[i].second;
            const int *conn = (*var.connectivity)[e];

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
                (var.mat->rho(e) + param.bc.wrinkler_delta_rho) * param.control.gravity * dz;

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
}


