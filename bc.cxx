#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "bc.hpp"


namespace {

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


bool is_on_boundary(const Variables &var, int node)
{
    uint flag = (*var.bcflag)[node];
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

        uint flag = (*var.bcflag)[i];
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
    if (param.bc.has_water_loading && param.control.gravity != 0) {

        const int top_bdry = bdry_order.find(BOUNDZ1)->second;
        const auto& top = var.bfacets[top_bdry];
        const auto& coord = *var.coord;
        // loops over all top facets
        for (std::size_t i=0; i<top.size(); ++i) {
            // this facet belongs to element e
            int e = top[i].first;
            // this facet is the f-th facet of e
            int f = top[i].second;
            const int *conn = (*var.connectivity)[e];

            // the outward-normal vector
            double normal[NDIMS];
            // the z-coordinate of the facet center
            double zcenter;

            normal_vector_of_facet(f, conn, *var.coord, normal, zcenter);

            const double sea_level = 0;
            const double sea_water_density = 1030;

            // below sea level?
            if (zcenter < sea_level) {
                double dz = sea_level - zcenter;
                double p = sea_water_density * param.control.gravity * dz;

                for (int j=0; j<NODES_PER_FACET; ++j) {
                    int n = conn[NODE_OF_FACET[f][j]];
                    for (int i=0; i<NDIMS; ++i) {
                        force[n][i] -= p * normal[i] / NODES_PER_FACET;
                    }
                }
            }
        }

    }

    // Wrinkler foundation for the bottom boundary
    if (param.bc.has_wrinkler_foundation && param.control.gravity != 0) {
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

            // the outward-normal vector
            double normal[NDIMS];
            // the z-coordinate of the facet center
            double zcenter;

            normal_vector_of_facet(f, conn, *var.coord, normal, zcenter);

            double dz = zcenter - (-param.mesh.zlength);
            double p = var.compensation_pressure -
                (var.mat->rho(e) + param.bc.wrinkler_delta_rho) * param.control.gravity * dz;

            // bottom support - Archimed force (normal to the surface)
            for (int j=0; j<NODES_PER_FACET; ++j) {
                int n = conn[NODE_OF_FACET[f][j]];
                for (int i=0; i<NDIMS; ++i) {
                    force[n][i] -= p * normal[i] / NODES_PER_FACET;
                }
            }
        }
    }
}


namespace {

    void simple_diffusion(const Variables& var, array_t& coord,
                          double surface_diffusivity)
    {
        /* Diffusing surface topography to simulate the effect of erosion and
         * sedimentation.
         */

        const int top_bdry = bdry_order.find(BOUNDZ1)->second;
        const auto& top = var.bfacets[top_bdry];

        const int_vec& top_nodes = var.bnodes[top_bdry];
        const std::size_t ntop = top_nodes.size();
        double_vec total_dx(var.nnode, 0);
        double_vec total_slope(var.nnode, 0);

        // loops over all top facets
        for (std::size_t i=0; i<top.size(); ++i) {
            // this facet belongs to element e
            int e = top[i].first;
            // this facet is the f-th facet of e
            int f = top[i].second;

            const int *conn = (*var.connectivity)[e];
            int n0 = (*var.connectivity)[e][NODE_OF_FACET[f][0]];
            int n1 = (*var.connectivity)[e][NODE_OF_FACET[f][1]];

#ifdef THREED
            int n2 = (*var.connectivity)[e][NODE_OF_FACET[f][2]];

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

            total_dx[n0] += projected_area;
            total_dx[n1] += projected_area;
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

                total_slope[n[j]] += slope * projected_area;
            }

            // std::cout << i << ' ' << n0 << ' ' << n1 << ' ' << n2 << "  "
            //           << projected_area << "  " << slope << '\n';
#else
            /* The 1D diffusion operation is implemented ad hoc, not using FEM
             * formulation (e.g. computing shape function derivation on the edges).
             */

            double dx = std::fabs(coord[n1][0] - coord[n0][0]);
            total_dx[n0] += dx;
            total_dx[n1] += dx;

            double slope = (coord[n1][1] - coord[n0][1]) / dx;
            total_slope[n0] -= slope;
            total_slope[n1] += slope;

            // std::cout << i << ' ' << n0 << ' ' << n1 << "  " << dx << "  " << slope << '\n';
#endif
        }

        double max_dh = 0;
        for (std::size_t i=0; i<ntop; ++i) {
            // we don't treat edge nodes specially, i.e. reflecting bc is used for erosion.
            int n = top_nodes[i];
            double dh = surface_diffusivity * var.dt * total_slope[n] / total_dx[n];
            coord[n][NDIMS-1] -= dh;
            max_dh = std::max(max_dh, std::fabs(dh));
            // std::cout << n << "  dh:  " << dh << '\n';
        }

        // std::cout << "max erosion / sedimentation rate (cm/yr):  "
        //           << max_dh / var.dt * 100 * YEAR2SEC << '\n';
    }


    void custom_surface_processes(const Variables& var, array_t& coord)
    {
        const int top_bdry = bdry_order.find(BOUNDZ1)->second;
        const int_vec& top_nodes = var.bnodes[top_bdry];
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

}


void surface_processes(const Param& param, const Variables& var, array_t& coord)
{
    switch (param.control.surface_process_option) {
    case 0:
        // no surface process
        break;
    case 1:
        simple_diffusion(var, coord, param.control.surface_diffusivity);
        break;
    case 101:
        custom_surface_processes(var, coord);
        break;
    default:
        std::cout << "Error: unknown surface process option: " << param.control.surface_process_option << '\n';
        std::exit(1);
    }
}


