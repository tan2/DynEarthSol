#include <algorithm>
#include <cstring>
#include <functional>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
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

    // std::cout << "tiny elements: ";
    // print(std::cout, tiny_elems);
    // std::cout << '\n';

    return tiny_elems.size();
}


void find_points_to_delete(const array_t &coord, const conn_t &connectivity,
                           const double_vec &volume, const int_vec &tiny_elems,
                           const array_t &old_coord, const int_vec &old_bcflag,
                           int_vec &to_delete)
{
    // collecting the nodes of tiny_elems
    int tiny_nelem = tiny_elems.size();
    array_t tiny_coord(tiny_nelem * NODES_PER_ELEM);
    conn_t tiny_conn(tiny_nelem);
    double_vec tiny_vol(tiny_nelem);
    int ii = 0;
    for (int ee=0; ee<tiny_nelem; ++ee) {
        int e = tiny_elems[ee];

        tiny_vol[ee] = volume[e];

        const int *conn = connectivity[e];
        for (int j=0; j<NODES_PER_ELEM; ++j) {
            int n = conn[j];
            tiny_conn[ee][j] = ii;

            for (int d=0; d<NDIMS; ++d) {
                tiny_coord[ii][d] = coord[n][d];
            }
            ii ++;
        }
    }

    Barycentric_transformation bary(tiny_coord, tiny_conn, tiny_vol);

    // find old nodes that are connected to tiny elements and are not on the boundary
    // (most of the nodes of tiny elements are newly inserted by the remeshing library)
    const int flag = BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | BOUNDZ0 | BOUNDZ1;
    for (int i=0; i<old_coord.size(); ++i) {
        // cannot delete boundary nodes
        if (old_bcflag[i] & flag) continue;

        const double *p = old_coord[i];
        for (int ee=0; ee<tiny_nelem; ++ee) {
            if (bary.is_inside_elem(p, ee)) {
                to_delete.push_back(i);
                break;
            }
        }
    }

    // std::cout << "old points to delete:\n";
    // print(std::cout, to_delete);
    // std::cout << '\n';
}


void delete_points(const int_vec &to_delete, const array_t &old_coord,
                   const segment_t &old_segment, double *points, int *segment)
{
    int nseg = old_segment.size();
    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = old_coord.size() - 1;
    // delete points from the end
    for (auto i=to_delete.rbegin(); i<to_delete.rend(); ++i) {
        // when a point is deleted, replace it with the last point
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        std::replace(segment, endsegment, end, *i);
        // std::cout << *i << " <- " << end << "\n";

        end --;
    }
}


void new_mesh(const Param &param, Variables &var,
              const array_t &old_coord, const conn_t &old_connectivity,
              const segment_t &old_segment, const segflag_t &old_segflag)
{
    int old_nnode = old_coord.size();
    int old_nseg = old_segment.size();

    // We don't want to refine large elements during remeshing,
    // so using the domain size as the max area
    double max_elem_size;
#ifdef THREED
    max_elem_size = param.mesh.xlength * param.mesh.ylength * param.mesh.zlength;
#else
    max_elem_size = param.mesh.xlength * param.mesh.zlength;
#endif
    const double vertex_per_polygon = 3;

    // new mesh
    int new_nnode, new_nelem, new_nseg;
    double *pcoord;
    int *pconnectivity, *psegment, *psegflag;
    points_to_new_mesh(param, old_nnode, old_coord.data(),
                       old_nseg, old_segment.data(), old_segflag.data(),
                       max_elem_size, vertex_per_polygon,
                       new_nnode, new_nelem, new_nseg,
                       pcoord, pconnectivity, psegment, psegflag);

    array_t new_coord(pcoord, new_nnode);
    conn_t new_connectivity(pconnectivity, new_nelem);

    // deleting (non-boundary) nodes to avoid having tiny elements
    double_vec new_volume(new_nelem);
    compute_volume(new_coord, new_connectivity, new_volume);

    int_vec tiny_elems;
    if (has_tiny_element(param, new_volume, tiny_elems)) {
        delete [] psegment;
        delete [] psegflag;

        int_vec to_delete;
        find_points_to_delete(new_coord, new_connectivity, new_volume,
                              tiny_elems, old_coord, *var.bcflag, to_delete);

        int q_nnode = old_nnode - to_delete.size();
        // create a copy of old_coord and old_segment
        double *qcoord = new double[old_coord.num_elements()];
        std::memcpy(qcoord, old_coord.data(), sizeof(double)*old_coord.num_elements());
        int *qsegment = new int[old_segment.num_elements()];
        std::memcpy(qsegment, old_segment.data(), sizeof(int)*old_segment.num_elements());

        delete_points(to_delete, old_coord, old_segment,
                      qcoord, qsegment);

        points_to_new_mesh(param, q_nnode, qcoord,
                           old_nseg, qsegment, old_segflag.data(),
                           max_elem_size, vertex_per_polygon,
                           new_nnode, new_nelem, new_nseg,
                           pcoord, pconnectivity, psegment, psegflag);

        delete [] qcoord;
        delete [] qsegment;
        new_coord.reset(pcoord, new_nnode);
        new_connectivity.reset(pconnectivity, new_nelem);
    }

    var.nnode = new_nnode;
    var.nelem = new_nelem;
    var.nseg = new_nseg;
    var.coord->steal_ref(new_coord);
    var.connectivity->steal_ref(new_connectivity);
    var.segment->reset(psegment, var.nseg);
    var.segflag->reset(psegflag, var.nseg);
}

} // anonymous namespace


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

    {
        // creating a "copy" of mesh pointer so that they are not deleted
        array_t old_coord;
        conn_t old_connectivity;
        segment_t old_segment;
        segflag_t old_segflag;
        old_coord.steal_ref(*var.coord);
        old_connectivity.steal_ref(*var.connectivity);
        old_segment.steal_ref(*var.segment);
        old_segflag.steal_ref(*var.segflag);

        new_mesh(param, var, old_coord, old_connectivity, old_segment, old_segflag);

        // interpolating fields
        nearest_neighbor_interpolation(var, old_coord, old_connectivity);
        barycentric_node_interpolation(var, old_coord, old_connectivity);

        // old_coord et al. are destroyed before exiting this block
    }

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

    std::cout << "  Remeshing finished.\n";
}


