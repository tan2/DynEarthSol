#include "iostream"

#include "constants.hpp"
#include "parameters.hpp"

#include "brc-interpolation.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "nn-interpolation.hpp"
#include "utils.hpp"
#include "remeshing.hpp"


bool bad_mesh_quality(const Param &param, const Variables &var)
{
    int worst_elem;
    double q = worst_elem_quality(*var.coord, *var.connectivity,
                                  *var.volume, worst_elem);
#ifdef THREED
    // normalizing q so that its magnitude is about the same in 2D and 3D
    q = std::pow(q, 1.0/3);
#endif
    std::cout << "Worst mesh quality = " << q << " at element #" << worst_elem << ".\n";
    if (q < param.mesh.min_quality) {
        return 1;
    }
    return 0;
}


void remesh(const Param &param, Variables &var)
{
    std::cout << "  Remeshing starts...\n";

    // setting up barycentric transformation

    // saving the mesh to .poly file

    // creating a "copy" of mesh pointer so that they are not deleted
    array_t old_coord;
    conn_t old_connectivity;
    segment_t old_segment;
    segflag_t old_segflag;
    old_coord.steal_ref(*var.coord);
    old_connectivity.steal_ref(*var.connectivity);
    old_segment.steal_ref(*var.segment);
    old_segflag.steal_ref(*var.segflag);

    delete var.coord;
    delete var.connectivity;
    delete var.segment;
    delete var.segflag;

    // XXX: modifying boundary if necessary

    // deleting (non-boundary) nodes to avoid small elements

    // saving the modified mesh to .poly file

    // new mesh
    int npoints = var.nnode;
    double *points = old_coord.data();
    int n_init_segments = var.nseg;
    int *init_segments = old_segment.data();
    int *init_segflags = old_segflag.data();

    // We don't want to refine large elements during remeshing,
    // so using the domain size as the max area
    double max_elem_size;
#ifdef THREED
    max_elem_size = param.mesh.xlength * param.mesh.ylength * param.mesh.zlength;
#else
    max_elem_size = param.mesh.xlength * param.mesh.zlength;
#endif
    double vertex_per_polygon = 3;
    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags,
                   max_elem_size, vertex_per_polygon);

    // interpolating fields
    nearest_neighbor_interpolation(var, old_coord, old_connectivity);
    barycentric_node_interpolation(var, old_coord, old_connectivity);

    // memory for new fields
    reallocate_variables(param, var);

    // updating other arrays
    delete var.bcflag;
    create_boundary_flags(var);
    for (int i=0; i<6; ++i) {
        var.bfacets[i].clear();
    }
    create_boundary_facets(var);
    delete var.support;
    create_support(var);
    delete var.egroups;
    create_elem_groups(var);

    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    free(old_coord.data());
    free(old_connectivity.data());
    free(old_segment.data());
    free(old_segflag.data());

    std::cout << "  Remeshing finished.\n";
}


