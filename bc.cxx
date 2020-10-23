#include <iostream>
#include <unordered_map>
#include <iomanip>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "markerset.hpp"


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
    return flag & BOUND_ANY;
}


double find_max_vbc(const BC &bc)
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


void create_boundary_normals(const Variables &var, double bnormals[nbdrytypes][NDIMS],
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
        if (var.bfacets[i].size() == 0) continue;

        for (auto j=var.bfacets[i].begin(); j<var.bfacets[i].end(); ++j) {
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

            if (j == var.bfacets[i].begin()) {
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
                    std::exit(1);
                }
            }
        }
    }

    for (int i=0; i<nbdrytypes; i++) {
        if (var.bfacets[i].size() == 0) continue;

        const double eps = 1e-15;
        for (int j=i+1; j<nbdrytypes; j++) {
            if (var.bfacets[j].size() == 0) continue;
            double *s = new double[NDIMS];  // intersection of two boundaries
                                            // whole-application lifetime, no need to delete manually
#ifdef THREED
            // quick path: both walls are vertical
            if (std::abs(var.bnormals[i][NDIMS-1]) < eps &&
                std::abs(var.bnormals[j][NDIMS-1]) < eps) {
                s[0] = s[1] = 0;
                s[NDIMS-1] = 1;
            }
            else {
                // cross product of 2 normal vectors
                s[0] = var.bnormals[i][1]*var.bnormals[j][2] - var.bnormals[i][2]*var.bnormals[j][1];
                s[1] = var.bnormals[i][2]*var.bnormals[j][0] - var.bnormals[i][0]*var.bnormals[j][2];
                s[2] = var.bnormals[i][0]*var.bnormals[j][1] - var.bnormals[i][1]*var.bnormals[j][0];
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

        //
        // X
        //
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
            case 7:
                v[0] = bc.vbc_val_x0;
                v[1] = 0;
                break;
#endif
            }
        }
        if (flag & BOUNDX1) {
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
            case 7:
                v[0] = bc.vbc_val_x1;
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
                v[0] = 0;
                v[1] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc.vbc_val_y0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc.vbc_val_y0;
                break;
            }
        }
        if (flag & BOUNDY1) {
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
                v[0] = bc.vbc_val_y1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y1;
                v[2] = 0;
                break;
            case 5:
                v[0] = bc.vbc_val_y1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 7:
                v[0] = 0;
                v[1] = bc.vbc_val_y1;
                break;
            }
        }
#endif

        //
        // N
        //
        for (int ib=iboundn0; ib<=iboundn3; ib++) {
            const double eps = 1e-15;
            const double *n = var.bnormals[ib]; // unit normal vector

            if (flag & (1 << ib)) {
                switch (var.vbc_types[ib]) {
                case 1:
                    if (flag == (1U << ib)) {  // ordinary boundary
                        double vn = 0;
                        for (int d=0; d<NDIMS; d++)
                            vn += v[d] * n[d];  // normal velocity

                        for (int d=0; d<NDIMS; d++)
                            v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                    }
                    else {  // intersection with another boundary
                        for (int ic=iboundx0; ic<ib; ic++) {
                            if (flag & (1 << ic)) {
                                if (var.vbc_types[ic] == 0) {
                                    double vn = 0;
                                    for (int d=0; d<NDIMS; d++)
                                        vn += v[d] * n[d];  // normal velocity

                                    for (int d=0; d<NDIMS; d++)
                                        v[d] += (var.vbc_values[ib] - vn) * n[d];  // setting normal velocity
                                }
                                else if (var.vbc_types[ic] == 1) {
                                    auto edge = var.edge_vectors.at(std::make_pair(ic, ib));
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
                        v[d] = var.vbc_values[ib] * n[d];  // v must be normal to n
                    break;
                }
            }
        }

        //
        // Z, must be dealt last
        //

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
            }
        }
        if (flag & BOUNDZ1) {
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
            }
        }
    }
}


void apply_stress_bcs(const Param& param, const Variables& var, array_t& force)
{
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

        const auto& bdry = var.bfacets[i];
        const auto& coord = *var.coord;
        // loops over all bdry facets
        for (int n=0; n<static_cast<int>(bdry.size()); ++n) {
            // this facet belongs to element e
            int e = bdry[n].first;
            // this facet is the f-th facet of e
            int f = bdry[n].second;
            const int *conn = (*var.connectivity)[e];

            // the outward-normal vector
            double normal[NDIMS];
            // the z-coordinate of the facet center
            double zcenter;

            normal_vector_of_facet(f, conn, *var.coord, normal, zcenter);

            double p;
            if (i==iboundz0 && param.bc.has_winkler_foundation) {
                // Winkler foundation for the bottom boundary
                p = var.compensation_pressure -
                    (var.mat->rho(e) + param.bc.winkler_delta_rho) *
                    param.control.gravity * (zcenter + param.mesh.zlength);
            }
            else if (i==iboundz1 && param.bc.has_water_loading) {
                // hydrostatic water loading for the surface boundary
                p = 0;
                if (zcenter < param.control.surf_base_level) {
                    // below sea level
                    const double sea_water_density = 1030;
                    p = sea_water_density * param.control.gravity * (param.control.surf_base_level - zcenter);
                }
            }
            else {
                // sidewalls
                p = ref_pressure(param, zcenter);
            }

            // lithostatc support - Archimed force (normal to the surface)
            for (int j=0; j<NODES_PER_FACET; ++j) {
                int n = conn[NODE_OF_FACET[f][j]];
                for (int d=0; d<NDIMS; ++d) {
                    force[n][d] -= p * normal[d] / NODES_PER_FACET;
                }
            }
        }
    }

    if (param.bc.has_elastic_foundation) {
        /* A restoration force on the bottom nodes proportional to total vertical displacement */
        for (auto i=var.bnodes[iboundz0].begin(); i<var.bnodes[iboundz0].end(); ++i) {
            int n = *i;
            force[n][NDIMS-1] -= param.bc.elastic_foundation_constant * ((*var.coord)[n][NDIMS-1] - (*var.coord0)[n][NDIMS-1]);
        }
    }
}


namespace {

    void simple_diffusion(const Variables& var, double_vec& dh)
    {
        /* Diffusing surface topography to simulate the effect of erosion and
         * sedimentation.
         */

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;

        const int top_bdry = iboundz1;
        const auto& top = var.bfacets[top_bdry];

        const std::size_t ntop = top_nodes.size();
        double_vec total_dx(var.nnode, 0);
        double_vec total_slope(var.nnode, 0);

        // loops over all top facets
#ifdef THREED
        for (std::size_t i=0; i<top.size(); ++i) {
            // this facet belongs to element e
            int e = top[i].first;
            // this facet is the f-th facet of e
            int f = top[i].second;

            const int *conn = (*var.connectivity)[e];
            int n0 = (*var.connectivity)[e][NODE_OF_FACET[f][0]];
            int n1 = (*var.connectivity)[e][NODE_OF_FACET[f][1]];

//#ifdef THREED
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
        for (std::size_t i=0; i<ntop-1;i++) {
            int n0 = top_nodes[i];
            int n1 = top_nodes[i+1];

            double dx = std::fabs(coord[n1][0] - coord[n0][0]);
            total_dx[n0] += dx;
            total_dx[n1] += dx;

            double slope = (coord[n1][1] - coord[n0][1]) / dx;
            total_slope[n0] -= slope;
            total_slope[n1] += slope;

            // std::cout << i << ' ' << n0 << ' ' << n1 << "  " << dx << "  " << slope << '\n';
#endif
        }

        for (std::size_t i=0; i<ntop; ++i) {
            // we don't treat edge nodes specially, i.e. reflecting bc is used for erosion.
            int n = top_nodes[i];
            double conv =  surfinfo.surf_diff * var.dt * total_slope[n] / total_dx[n];
            if ( coord[n][1] >  surfinfo.base_level && conv > 0.) {
                dh[i] -= surfinfo.diff_ratio_terrig * conv;
            } else if ( coord[n][1] < surfinfo.base_level && conv < 0. ) {
                dh[i] -= surfinfo.diff_ratio_marine * conv;
            } else{
                dh[i] -= conv;
            }
        }
    }


    void custom_surface_processes(const Variables& var, array_t& coord) {
        const int top_bdry = iboundz1;
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

    void get_surface_info(const Variables& var, \
        double_vec& top_base, double_vec& top_slope, double_vec& top_deg, \
        double_vec& top_depth) {

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

            double slope = (coord[n1][1] - coord[n0][1]) / dx;
            top_slope[i] += slope;
            top_slope[i+1] += slope;
        }

        for (std::size_t i=0; i<ntop; i++) {
            top_depth[i] = surfinfo.base_level - coord[top_nodes[i]][1];
            top_deg[i] = atan(top_slope[i]);
            top_base[i] *= 0.5;
        }

/*
        for (std::size_t i=0; i<ntop-1; i++) {
            top_slope[i] = ( coord[top_nodes[i+1]][1] - coord[top_nodes[i]][1] );
            top_base[i] = ( coord[top_nodes[i+1]][0] - coord[top_nodes[i]][0] );
            top_slope[i] /= top_base[i];
            top_deg[i] = atan(top_slope[i]);
        }


        // calculate base of node
        top_base[ntop-1] = top_base[ntop-2] / 2.;
        for (std::size_t i=ntop-2; i>0; i--)
            top_base[i] = ( top_base[i] + top_base[i-1] ) / 2.;
        top_base[0] = top_base[0] / 2.;
*/
    }

    void get_basin_info(const Variables& var, double_vec& top_depth, \
        std::vector<bool>& if_source, int_vec2D& source_pair) {

        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();

        int_vec if_land(ntop,0);
        int_vec imouth_loc;
        int_vec left_mouse, right_mouse;
        int imouth_left, imouth_right;
        int_vec mouths(ntop,0);

        // find land:1 and sea:0
        for (std::size_t i=0; i<ntop; i++)
            if (top_depth[i] <= 0.) if_land[i] = 1;

        // find river mouth
        for (std::size_t i=0; i<ntop-1; i++)
            if (if_land[i] + if_land[i+1] == 1) {
                if (if_land[i] == 1) {
                    imouth_loc.push_back(i);
                    left_mouse.push_back(i);
                    mouths[i] += 1;
                }
                else {
                    imouth_loc.push_back(i+1);
                    right_mouse.push_back(i+1);
                    mouths[i+1] += 2;
                }
            }
        
        // find basin sets
        int left_mouth_num = left_mouse.size();
        int right_mouth_num = right_mouse.size();

        if ( var.steps%10000 == 0 ) {
            for (std::size_t i=0;i<ntop;i++)
                printf("%d",if_land[i]);
            printf("\n");
        }

        int left_source_point;
        int right_source_point;

        if (left_mouth_num > 0) {
            if_source[0] = true;
            left_source_point = left_mouse[0];
        }
        if (right_mouth_num > 0) {
            if_source[1] = true;
            right_source_point = right_mouse[right_mouth_num-1];
        }

        // find basin pairs
        if (if_source[0] == true) {
            if (right_mouth_num != 0 and left_source_point < right_mouse[0]) {
                source_pair[0][0] = left_source_point;
                source_pair[0][1] = right_mouse[0];
            }
            else {
                source_pair[0][0] = left_source_point;
                source_pair[0][1] = ntop-3;
            }
        }

        if (if_source[1] == true ) {
            if (left_mouth_num != 0 and left_mouse[left_mouth_num-1] < right_source_point) {
                source_pair[1][0] = left_mouse[left_mouth_num-1];
                source_pair[1][1] = right_source_point;
            }
            else {
                source_pair[1][0] = 2;
                source_pair[1][1] = right_source_point;                
            }
        }
    }


    void simple_deposition(const Param& param,const Variables& var, double_vec& dh) {

#ifdef THREED
        // not ready for 3D
        std::cout << "3D deposition of sediment processes is not ready yet.";
        exit(168);
//        std::cout << "Press enter to continue ...";
//        std::cin.get();
#endif
        const array_t& coord = *var.coord;
        const SurfaceInfo& surfinfo = var.surfinfo;
        const int_vec& top_nodes = *surfinfo.top_nodes;
        const std::size_t ntop = top_nodes.size();

        double_vec top_base(ntop,0.);
        double_vec top_slope(ntop,0.);
        double_vec top_deg(ntop,0.);
        double_vec top_depth(ntop,0.);

        std::vector<bool> if_source(2,false);
        double_vec dh_terrig(ntop,0.);
        double_vec dhacc_tmp(ntop,0.);

        if ( var.steps%10000 == 0 )
            printf("** Sedimentation report:\n");
        get_surface_info(var,top_base,top_slope,top_deg,top_depth);

//******************************************************************
//              sedimentation by terrigenous source
//******************************************************************

        // The Pearl River Shibao et al. (2007)
        // 5.e7 ton/yr ~= 1. m^3/s
        // assuming the width is 50 km --> 2.e-5 m^2/s
        double max_sedi_vol = param.control.surf_src_vol;
        if (max_sedi_vol != 1.)
            max_sedi_vol = max_sedi_vol / 50.e3 * var.dt;
        else
            max_sedi_vol = param.control.surf_src_area * var.dt;

        double vol_ratio = 10;
        int ntry = 200;
        int_vec2D source_pair(2,int_vec(2,0));
        double basin_vol, dh_tmp, mul_conv;

        get_basin_info(var,top_depth,if_source, source_pair);

        int_vec start(2,0);
        int_vec end(2,0);
        int_vec sign(2,0);
        int_vec width(2,0);
        int_vec isedi(2,0);

        int_vec2D loccation_index(2,int_vec(ntop,0));
        double_vec sedi_vol(2,0.);
        double_vec unit_vol(2,0.);
        std::vector<bool> if_space_limited(2,false);
        std::vector<bool> if_slope_limited(2,false);

        start[0] = source_pair[0][0];
        end[0] = source_pair[0][1];
        sign[0] = 1;

        start[1] = source_pair[1][1];
        end[1] = source_pair[1][0]; 
        sign[1] = -1;

        // deal with the sedimentation for the aspect of coastal source
        for (int i=0; i<2; i++) {

            if (! if_source[i]) continue;
            // the source is too close to boundary
            if (abs(start[i]-int(ntop)/2) >= int(ntop)/2) {
                if_source[i] = false;
                continue;
            }
            width[i] = abs(end[i] - start[i]);

            loccation_index[i][0] = start[i];
            for (int j=1; j<width[i];j++)
                loccation_index[i][j] = loccation_index[i][j-1] + sign[i];

            basin_vol = 0;
            for (int ind=0; ind<width[i]; ind++) {
                int j = loccation_index[i][ind];
                basin_vol += top_depth[j] * top_base[j];
            }
            // the space of basin is too small
//            if (basin_vol < max_sedi_vol / 10.) if_source[i] = false;
        }

        //
        for (int k=0; k<10;k++) {

            for (int i=0; i<2; i++) {
                if (! if_source[i]) continue;

                unit_vol[i] = std::min(max_sedi_vol / 10., max_sedi_vol - sedi_vol[i]);
                bool if_finish = false;
                int iloop = 0;
                double dist;
                double_vec depo_vol(2,0.);

                do {
                    iloop++;
                    
                    // calculate the dh of basin
                    for (int ind=1; ind<width[i]; ind++) {
                        int j = loccation_index[i][ind];
                        int pj = loccation_index[i][ind-1];
                        // stop propagate if slope is increasing
                        if ( top_slope[j - 1 + i ] * sign[i] > 0. ) {
                            if_slope_limited[i] = true;
                            break;
                        }

                        dist = abs(coord[top_nodes[j]][0] - coord[top_nodes[start[i]]][0]);

                        dh_tmp = 8.e-9 * pow(1.25,iloop) * (unit_vol[i] - depo_vol[i]) * dist;

                        if ( dh_tmp != dh_tmp ) {
                            printf("dh_tmp is NaN.\n");
                            printf("%e\t%e\t%f\n",dh_tmp, (unit_vol[i] - depo_vol[i]), dist);
                        }

                        dhacc_tmp[j] += dh_tmp;

                        if (dhacc_tmp[j] > top_depth[j])
                            dhacc_tmp[j] = top_depth[j] + 0.1;

                        depo_vol[i] += top_base[j] * dh_tmp;

                        if (depo_vol[i] >= unit_vol[i]) {
                            if_finish = true;
                            break;
                        }
                    }
                    

                    // move coast if is full.
                    bool if_shift = false;
                    if (dhacc_tmp[start[i]+sign[i]] >= top_depth[start[i]+sign[i]]) {
                        start[i] += sign[i];
                        if_shift = true;
                        if (start[i] == end[i]) {
                            if_space_limited[i] = true;
                            break;
                        }
                    }

                    // if the coast is shifted, change the coast index
                    if (if_shift) {
                        width[i] = abs(end[i] - start[i]);
                        loccation_index[i][0] = start[i];
                        for (int j=1; j<width[i]; j++)
                            loccation_index[i][j] = loccation_index[i][j-1] + sign[i];
                    }

                } while (unit_vol[i] - depo_vol[i] > 1.e-4 && !if_finish && iloop <= ntry);

                isedi[i] += iloop;
                sedi_vol[i] += depo_vol[i];
            }

            if (if_space_limited[0] || if_space_limited[1]) break;

            for (std::size_t j=0; j<ntop; j++) {
                dh_terrig[j] += dhacc_tmp[j];
                dhacc_tmp[j] = 0.;
            }
            if (sedi_vol[0] >= max_sedi_vol || sedi_vol[1] >= max_sedi_vol) break;


        }

        if ( var.steps%10000 == 0 ) {
            for (int i=0;i<2;i++) {
                if (if_space_limited[i]) printf("   Space limited at %d\n",i);
                if (if_slope_limited[i]) printf("   Slope limited at %d\n",i);
                if (!if_source[i]) {
                    printf("    No source from %d\n",i);
                } else {
                    printf("    Side %d: Location: %10.2f km (Index: %5d - %5d). Sediment: %10.2f m^2 (max: %10.2f) Loop: %5d\n",\
                        i,coord[top_nodes[start[i]]][0]/1000.,start[i],end[i],sedi_vol[i], max_sedi_vol,isedi[i]);
                }
            }

            if (if_space_limited[0] || if_space_limited[1])
                printf("    Space of basin is not enough for sediment . Do next round. %5d/%5d\n",isedi[0],isedi[1]);
        }

        for (std::size_t i=0; i<ntop; i++)
            dh[i] += dh_terrig[i];

//******************************************************************
//              sedimentation by suspended source
//******************************************************************

        double ddh = surfinfo.depo_universal * var.dt;
        for (std::size_t i=0; i<ntop; i++)
            // if below the base level
            if (top_depth[i] > 0.) dh[i] += ddh;
    }

}


//void surface_processes(const Param& param, const Variables& var, array_t& coord, double_vec& plstrain, \
//                      SurfaceInfo& surfinfo, std::vector<MarkerSet*> &markersets, int_vec2D& elemmarkers)
void surface_processes(const Param& param, const Variables& var, array_t& coord, \
                       SurfaceInfo& surfinfo, std::vector<MarkerSet*> &markersets, int_vec2D& elemmarkers)
{
    int ntop = surfinfo.top_nodes->size();
    double_vec dh(ntop,0.);
    
    switch (param.control.surface_process_option) {
    case 0:
        // no surface process
        break;
    case 1:
        simple_diffusion(var, dh);
        break;
    case 101:
        custom_surface_processes(var, coord);
        break;
    case 102:
        simple_diffusion(var, dh);
        simple_deposition(param, var, dh);
        break;
    default:
        std::cout << "Error: unknown surface process option: " << param.control.surface_process_option << '\n';
        std::exit(1);
    }

    if ( var.steps%10000 == 0 ) {
        double max_dh = 0., min_dh = 0.;
        for (int i=0;i<ntop;i++) {
            max_dh = std::max(max_dh, dh[i]);
            min_dh = std::min(min_dh, dh[i]);
            // std::cout << n << "  dh:  " << dh << '\n';
        }
        std::cout << "max erosion / sedimentation rate (mm/yr):  "
                    << min_dh / var.dt * 1000 * YEAR2SEC << " / "
                    << max_dh / var.dt * 1000 * YEAR2SEC << '\n';
    }

    // go through all surface nodes and abject all surface node by dh
    for (std::size_t i=0; i<ntop; i++) {
        // get global index of node
        int n = (*surfinfo.top_nodes)[i];
        // update coordinate via dh
        (coord)[n][NDIMS-1] += dh[i];
        // update dhacc for marker correction
        (*surfinfo.dhacc)[n] += dh[i];

        // go through connected elements
        for (std::size_t j=0; j<(*surfinfo.node_and_elems)[i].size(); j++) {
            // get local index of surface element
            int e = (*surfinfo.node_and_elems)[i][j];
            // get global index of element
            int eg = (*surfinfo.top_facet_elems)[e];
            // get local index of node in connected element
            int ind = (*surfinfo.arcelem_and_nodes_num)[e][i];
            // update edhacc of connected elements
            (*surfinfo.edhacc)[eg][ind] += dh[i];
        }
    }

    if (!(var.steps % param.mesh.quality_check_step_interval)) {
        // correct surface element value.
//        correct_surface_element(var,param.mat.mattype_sed,plstrain);
        // correct surface marker.
        markersets[0]->correct_surface_marker(var);
        std::fill(surfinfo.dhacc->begin(), surfinfo.dhacc->end(), 0.);
        // set marker of sediment.
        markersets[0]->set_sediment_marker(var,param.mat.mattype_sed,*surfinfo.edhacc, elemmarkers);

    }
}
