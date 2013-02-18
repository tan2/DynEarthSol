#include <algorithm>
#include <functional>
#include <iostream>

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

namespace {


bool has_tiny_element(const Param &param, const double_vec &volume,
                      int_vec &tiny_elems)
{
    const double smallest_vol = param.mesh.smallest_size * std::pow(param.mesh.resolution, NDIMS);

    for (int e=0; e<volume.size(); e++) {
        if (volume[e] < smallest_vol)
            tiny_elems.push_back(e);
    }

    // print(std::cout, tiny_elems);

    return tiny_elems.size();
}


void find_points_to_delete(const conn_t &connectivity, const int_vec &tiny_elems,
                           int_vec &to_delete)
{
    // find nodes that are connected to tiny elements and count the # of connection
    std::map<int, int> common_nodes;
    for (int i=0; i<tiny_elems.size(); ++i) {
        int e = tiny_elems[i];
        const int *conn = connectivity[e];
        for (int j=0; j<NODES_PER_ELEM; ++j) {
            int node = conn[j];
            common_nodes[node] += 1;
        }
    }

    // mark nodes that are connected to two (or more) tiny elements
    for (auto i=common_nodes.begin(); i!=common_nodes.end(); ++i) {
        if (i->second > 1)
            to_delete.push_back(i->first);
    }
}


void delete_points(const int_vec &bcflag, int_vec &to_delete,
                   int &npoints, double *points,
                   int &n_init_segments, int *init_segments, int *init_segflags)
{
    // will delete points in descending order
    std::sort(to_delete.begin(), to_delete.end(), std::greater<int>());

    std::cout << "points to delete:\n";
    print(std::cout, to_delete);
    std::cout << '\n';

    std::cout << "segment before:\n";
    print(std::cout, init_segments, n_init_segments*NODES_PER_FACET);
    std::cout << '\n';

    std::cout << npoints << ' ' << to_delete.size() << '\n';

    for (auto i=to_delete.begin(); i!=to_delete.end(); ++i) {
        std::cout << "deleting " << *i << " with flag " << bcflag[*i] << '\n';

        // XXX: cannot delete boundary nodes, otherwise the segments and segflags
        // will be messed up
        if (bcflag[*i]) continue;

        // when a point is deleted, replace it with the last point
        npoints--;
        std::cout << "  replaced by " << npoints << '\n';
        for (int d=0; d<NDIMS; ++d)
            points[(*i)*NDIMS + d] = points[npoints*NDIMS + d];

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        int *endseg = init_segments + n_init_segments * NODES_PER_FACET;
        std::replace(init_segments, endseg, npoints, *i);
    }

    std::cout << "segment after:\n";
    print(std::cout, init_segments, n_init_segments*NODES_PER_FACET);
    std::cout << '\n';

    std::cout << npoints << '\n';
}

}


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

    std::cout << "segment:\n";
    print(std::cout, *var.segment);
    std::cout << '\n';
    std::cout << "segflag:\n";
    print(std::cout, *var.segflag);
    std::cout << '\n';

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
    const double vertex_per_polygon = 3;

    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags,
                   max_elem_size, vertex_per_polygon);

    // deleting (non-boundary) nodes to avoid having tiny elements
    double_vec new_volume(var.nelem);
    compute_volume(*var.coord, *var.connectivity, new_volume);

    int_vec to_delete;
    int_vec tiny_elems;
    if (has_tiny_element(param, new_volume, tiny_elems)) {
        delete var.bcflag;
        create_boundary_flags(var);

        find_points_to_delete(*var.connectivity, tiny_elems, to_delete);
        delete_points(*var.bcflag, to_delete, var.nnode, var.coord->data(),
                      var.nseg, var.segment->data(), var.segflag->data());

        points_to_mesh(param, var, npoints, points,
                       n_init_segments, init_segments, init_segflags,
                       max_elem_size, vertex_per_polygon);
    }

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

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    // TODO: using edvoldt and volume to get volume_old
    std::copy(var.volume->begin(), var.volume->end(), var.volume_old->begin());
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.volume_n, *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    // the following variables need to be re-computed only because we are
    // outputing right after remeshing
    update_strain_rate(var, *var.strain_rate);
    update_force(param, var, *var.force);

    free(old_coord.data());
    free(old_connectivity.data());
    free(old_segment.data());
    free(old_segflag.data());

    std::cout << "  Remeshing finished.\n";
}


