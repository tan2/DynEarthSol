#include <algorithm>
#include <cstring>
#include <functional>
#include <iostream>
#include <set>

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
#include "markerset.hpp"
#include "remeshing.hpp"

namespace {


// equilateral triangle area = 0.433*s^2, equilateral tetrahedron volume = 0.118*s^3
#ifdef THREED
const double sizefactor = 0.118;
#else
const double sizefactor = 0.433;
#endif

const int DELETED_FACET = -1;
const int DEBUG = 0;

bool is_boundary(uint flag)
{
    return flag & (BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | BOUNDZ0 | BOUNDZ1);
}


bool is_bottom(uint flag)
{
    return flag & BOUNDZ0;
}


bool is_corner(uint flag)
{
    const std::set<uint> corner_set =
        {
#ifdef THREED
            BOUNDX0 | BOUNDY0 | BOUNDZ0,
            BOUNDX0 | BOUNDY0 | BOUNDZ1,
            BOUNDX0 | BOUNDY1 | BOUNDZ0,
            BOUNDX0 | BOUNDY1 | BOUNDZ1,
            BOUNDX1 | BOUNDY0 | BOUNDZ0,
            BOUNDX1 | BOUNDY0 | BOUNDZ1,
            BOUNDX1 | BOUNDY1 | BOUNDZ0,
            BOUNDX1 | BOUNDY1 | BOUNDZ1,
#else
            BOUNDX0 | BOUNDZ0,
            BOUNDX0 | BOUNDZ1,
            BOUNDX1 | BOUNDZ0,
            BOUNDX1 | BOUNDZ1,
#endif
        };

    uint f = flag & (BOUNDX0 | BOUNDX1 |
                     BOUNDY0 | BOUNDY1 |
                     BOUNDZ0 | BOUNDZ1);

    if (!f) return 0;

    if (corner_set.find(f) == corner_set.end()) return 0;
    return 1;
}


bool is_bottom_corner(uint flag)
{
    if ((flag & BOUNDZ0) && is_corner(flag)) return 1;
    return 0;
}


void flatten_bottom(const uint_vec &old_bcflag, double *qcoord,
                    double bottom, int_vec &points_to_delete, double min_dist)
{
    // find old nodes that are on or close to the bottom boundary

    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_bottom(flag)) {
            // restore edge nodes to initial depth
            qcoord[i*NDIMS + NDIMS-1] = bottom;
        }
        else if (! is_boundary(flag) &&
                 std::fabs(qcoord[i*NDIMS + NDIMS-1] - bottom) < min_dist) {
            points_to_delete.push_back(i);
        }
    }
}


void new_bottom(const uint_vec &old_bcflag, double *qcoord,
                double bottom_depth, int_vec &points_to_delete, double min_dist,
                int *segment, int *segflag, int nseg)
{
    /* deleting nodes that are on or close to the bottom boundary,
     * excluding nodes on the side walls
     */

    int_vec bottom_corners;
    for (std::size_t i=0; i<old_bcflag.size(); ++i) {
        uint flag = old_bcflag[i];
        if (is_bottom(flag)) {
            if(is_bottom_corner(flag))
                bottom_corners.push_back(i);
            else
                points_to_delete.push_back(i);
        }
        else if (! is_boundary(flag) &&
                 std::fabs(qcoord[i*NDIMS + NDIMS-1] - bottom_depth) < min_dist) {
            points_to_delete.push_back(i);
        }
    }

    if (DEBUG) {
        std::cout << "bottom points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
        std::cout << "segment before delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag before delete: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }

    // must have 2 corners in 2D, 4 corners in 3D
    if (bottom_corners.size() != (2 << (NDIMS-2))) {
        std::cerr << "Error: cannot find all bottom corners before remeshing. n_bottom_corners = "
                  << bottom_corners.size() << '\n';
        std::cout << "bottom corners: ";
        print(std::cout, bottom_corners);
        std::cout << '\n';
        std::exit(1);
    }

    // move the corners to the same depth
    for (std::size_t i=0; i<bottom_corners.size(); i++) {
        int n = bottom_corners[i];
        qcoord[n*NDIMS + NDIMS-1] = bottom_depth;
    }

    // mark all bottom facets nodes as deleted
    for (int i=0; i<nseg; ++i) {
        if (static_cast<uint>(segflag[i]) == BOUNDZ0) {
            for (int j=0; j<NODES_PER_FACET; j++)
                segment[i*NODES_PER_FACET + j] = DELETED_FACET;
        }
    }

    // create new bottom facets from corner nodes
    // XXX: Assuming square box, 1 facet (segment) in 2D, 2 facets in 3D
    // XXX: Assuming the order of bottom nodes in 3D is
    //      0 -- 1
    //      |    |
    //      2 -- 3
    for (int i=0, nfacets=0, offset=0; i<nseg, nfacets<(NDIMS-1); ++i) {
        if (static_cast<uint>(segflag[i]) == BOUNDZ0) {
            for (int j=0; j<NODES_PER_FACET; j++)
                segment[i*NODES_PER_FACET + j] = bottom_corners[offset + j];

            segflag[i] = BOUNDZ0;
            nfacets ++;
            offset ++;
        }
    }

    if (DEBUG) {
        std::cout << "bottom corners: ";
        print(std::cout, bottom_corners);
        std::cout << '\n';
        std::cout << "segment with new bottom: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag with new bottom: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }
}


void assemble_facet_polygons(const Variables &var, const array_t &old_coord, int_vec (&facet_polygons)[6])
{
    /* facet_polygons[i] contains a list of vertex which enclose the i-th boundary facet */
#ifdef THREED
    // 3d remeshing might need additional info on edge connection
    const int nedges = 12;
    int_vec edgenodes[nedges];

    /* Expanded diagram of the brick, # vertex, E# edge, F# facet
     *         4----E3-----7
     *         |           |
     *         |    F3     |
     *         |           |
     *   4-E10-5----E1-----6-E11-7----E3-----4
     * E6|     |E4         |E5   |E7         |E6
     *   | F0  |    F4     | F1  |     F5    |
     *   |     |           |     |           |
     *   0-E8--1----E0-----2-E9--3----E2-----0
     *         |           |
     *         |    F3     |
     *         |           |
     *         0----E2-----3
     */


    // WARNING: this part has hardcoded numbers for boundary and edges
    // MIGHT BREAK when the domain is not a brick
    for (int j=0; j<2; ++j) {
        const int_vec &bnode = var.bnodes[j];
        for (std::size_t i=0; i<bnode.size(); ++i) {
            uint f = (*var.bcflag)[bnode[i]];
            if (f & BOUNDZ0) {
                edgenodes[j+4].push_back(bnode[i]);
            }
            if (f & BOUNDZ1) {
                edgenodes[j+6].push_back(bnode[i]);
            }
            if (f & BOUNDY0) {
                edgenodes[j+8].push_back(bnode[i]);
            }
            if (f & BOUNDY1) {
                edgenodes[j+10].push_back(bnode[i]);
            }
        }
    }
    for (int j=2; j<4; ++j) {
        const int_vec &bnode = var.bnodes[j];
        for (std::size_t i=0; i<bnode.size(); ++i) {
            uint f = (*var.bcflag)[bnode[i]];
            if (f & BOUNDZ0) {
                edgenodes[j-2].push_back(bnode[i]);
            }
            if (f & BOUNDZ1) {
                edgenodes[j+0].push_back(bnode[i]);
            }
        }
    }

    // local cmp functor
    struct cmp {
        const int d;
        const array_t &coord;
        cmp (const array_t &coord_, int dim) : coord(coord_), d(dim) {};
        int operator()(const int a, const int b) {return coord[a][d] < coord[b][d];}
    };

    // ordering the edge nodes by sorting the coordinate along the the corresponding dimension
    // edge 0~3 along X; edge 4~7 along Y; edge 8~11 along Z
    for (int i=0; i<nedges; ++i) {
        int_vec &enodes = edgenodes[i];
        int dim = i / 4;
        std::sort(enodes.begin(), enodes.end(), cmp(old_coord, dim));
    }

    if (DEBUG > 1) {
        for (int j=0; j<nedges; ++j) {
            std::cout << "edge " << j << '\n';
            print(std::cout, edgenodes[j]);
            std::cout << '\n';
        }
    }

    // polygons that enclosing each boundary facets (orientation is outward normal)
    const int edgelist[6][4] = {{ 8, 6,10, 4},
                                { 5,11, 7, 9},
                                { 0, 9, 2, 8},
                                {10, 3,11, 1},
                                { 4, 1, 5, 0},
                                { 2, 7, 3, 6}};

    for (int n=0; n<6; ++n) {
        int j;
        j = edgelist[n][0];
        // ending at end()-1 to avoid duplicating corner nodes
        for (auto i=edgenodes[j].begin(); i!=edgenodes[j].end()-1; ++i) facet_polygons[n].push_back(*i);
        j = edgelist[n][1];
        for (auto i=edgenodes[j].begin(); i!=edgenodes[j].end()-1; ++i) facet_polygons[n].push_back(*i);
        j = edgelist[n][2];
        // using reversed iterator to maintain orientation
        for (auto i=edgenodes[j].rbegin(); i!=edgenodes[j].rend()-1; ++i) facet_polygons[n].push_back(*i);
        j = edgelist[n][3];
        for (auto i=edgenodes[j].rbegin(); i!=edgenodes[j].rend()-1; ++i) facet_polygons[n].push_back(*i);
    }

    if (DEBUG > 1) {
        for (int j=0; j<6; ++j) {
            std::cout << "facet polygon nodes " << j << '\n';
            print(std::cout, facet_polygons[j]);
            std::cout << '\n';
        }
    }

#endif
    return;
}


void find_tiny_element(const Param &param, const double_vec &volume,
                       int_vec &tiny_elems)
{
    const double smallest_vol = param.mesh.smallest_size * sizefactor * std::pow(param.mesh.resolution, NDIMS);

    for (std::size_t e=0; e<volume.size(); e++) {
        if (volume[e] < smallest_vol)
            tiny_elems.push_back(e);
    }

    if (DEBUG) {
        std::cout << "tiny elements: ";
        print(std::cout, tiny_elems);
        std::cout << '\n';
    }
}


void find_points_of_tiny_elem(const array_t &coord, const conn_t &connectivity,
                              const double_vec &volume, const int_vec &tiny_elems,
                              int npoints, const double *old_points,
                              const uint_vec &old_bcflag, int_vec &points_to_delete,
                              bool excl_func(uint))
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

    // find old nodes that are connected to tiny elements and are not excluded
    // (most of the nodes of tiny elements are newly inserted by the remeshing library)
    for (int i=0; i<npoints; ++i) {
        // excluded nodes
        if (excl_func(old_bcflag[i])) continue;

        const double *p = old_points + i*NDIMS;
        for (int ee=0; ee<tiny_nelem; ++ee) {
            if (bary.is_inside_elem(p, ee)) {
                points_to_delete.push_back(i);
                break;
            }
        }
    }

    if (DEBUG) {
        std::cout << "points of tiny elements: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
    }
}


void delete_points(int_vec &points_to_delete, int &npoints,
                   int nseg, double *points, int *segment)
{
    /* sort points_to_delete and remove duplicate */
    std::sort(points_to_delete.begin(), points_to_delete.end());
    auto last = std::unique(points_to_delete.begin(), points_to_delete.end());
    points_to_delete.resize(last - points_to_delete.begin());

    if (DEBUG) {
        std::cout << "old points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
    }

    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = npoints - 1;

    // delete points from the end
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
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
    npoints -= points_to_delete.size();
}


void delete_points_and_merge_segments(const int_vec &points_to_delete, int &npoints,
                                      int nseg, double *points, int *segment,
                                      uint_vec &bcflag, double min_length)
{
#ifdef THREED
    std::cerr << "delete_points_and_merge_segments() doesn't work in 3D!\n";
    std::exit(1);
#endif

    int *endsegment = segment + nseg * NODES_PER_FACET;

    int end = npoints - 1;

    // delete points from the end
    for (auto i=points_to_delete.rbegin(); i<points_to_delete.rend(); ++i) {
        if (DEBUG) {
            std::cout << " deleting point " << *i << " replaced by point " << end << '\n';
        }

        // if the deleted point is a segment point, merge the neighboring two segments
        uint flag = bcflag[*i];
        if (is_boundary(flag)) {
            // segment and points:  aa------a = b-------bb

            int *a = std::find(segment, endsegment, *i);
            int *b = std::find(a+1, endsegment, *i);
            if (b == endsegment) {
                std::cerr << "Error: segment array is corrupted when merging segment!\n";
                std::exit(1);
            }

            // a could be the either first or second node of the segment,
            // which will affect the location of aa in the segment array.
            bool not_first;
            not_first = (a - segment) % NODES_PER_FACET;
            int *aa = not_first? (a - 1) : (a + 1);
            not_first = (b - segment) % NODES_PER_FACET;
            int *bb = not_first? (b - 1) : (b + 1);

            // if the length of the two segments are too long
            // do not delete this point
            double la2, lb2;
            la2 = dist2(points + (*a)*NDIMS, points + (*aa)*NDIMS);
            lb2 = dist2(points + (*b)*NDIMS, points + (*bb)*NDIMS);
            if (la2 > min_length*min_length && lb2 > min_length*min_length) {
                if (DEBUG) {
                    std::cout << " the segments of point " << *i << " have length^2 "
                              << la2 << ", " << lb2 << " -- skip deletion."<< '\n';
                }
                continue;
            }

            // merge the two segments
            *a = *bb;
            *bb = DELETED_FACET;
            *b = DELETED_FACET;

            if (DEBUG) {
                std::cout << "a: " << (a-segment) << "  b: " << (b-segment)
                          << " not_first? " << not_first << " bb: " << (bb-segment)
                          << " = " << *bb << '\n';
                std::cout << "segment after merging: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        // when a point is deleted, replace it with the last point
        flag = bcflag[end];
        bcflag[*i] = flag;
        for (int d=0; d<NDIMS; ++d) {
            points[(*i)*NDIMS + d] = points[end*NDIMS + d];
        }

        // if the last point is also a segment point, the segment point index
        // needs to be updated as well
        if (is_boundary(flag)) {
            std::replace(segment, endsegment, end, *i);
            // std::cout << *i << " <- " << end << "\n";
            if (DEBUG) {
                std::cout << "segment after replace: ";
                print(std::cout, segment, nseg*NODES_PER_FACET);
                std::cout << '\n';
            }
        }

        end --;
        npoints --;
    }

    /* PS: We don't check whether the merged segments are belong the same boundary or not,
     *     e.g., one segment is on the top boundary while the other segment is on the left
     *     boundary. The check is performed in delete_facets() instead.
     */
}


void delete_points_and_merge_facets()
{
#ifndef THREED
    std::cerr << "delete_points_and_merge_facets() doesn't work in 2D!\n";
    std::exit(1);
#endif



}


void delete_points_on_boundary(const int_vec &points_to_delete, int &npoints,
                               int nseg, double *points, int *segment,
                               uint_vec &bcflag, double min_length)
{
    if (DEBUG) {
        std::cout << "old points to delete: ";
        print(std::cout, points_to_delete);
        std::cout << '\n';
        std::cout << "segment before delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
    }

#ifdef THREED
    delete_points_and_merge_facets();
#else
    delete_points_and_merge_segments(points_to_delete, npoints, nseg,
                                     points, segment, bcflag, min_length);
#endif

    if (DEBUG) {
        std::cout << "segment after  delete: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
    }
}


void delete_facets(int &nseg, int *segment, int *segflag)
{
    // delete facets from the end
    for (int i=nseg-1; i>=0; i--) {
        if (segment[i*NODES_PER_FACET] == DELETED_FACET) {
            // safety check
            if (segment[i*NODES_PER_FACET + 1] != DELETED_FACET
#ifdef THREED
                || segment[i*NODES_PER_FACET + 2] != DELETED_FACET
#endif
                ) {
                std::cerr << "Error: segment array is corrupted before delete_facets()!\n";
                print(std::cerr, segment, nseg*NODES_PER_FACET);
                std::exit(1);
            }

            // replace deleted segment with the last segment
            for (int j=0; j<NODES_PER_FACET; ++j) {
                segment[i*NODES_PER_FACET + j] = segment[(nseg-1)*NODES_PER_FACET + j];
            }
            segflag[i] = segflag[nseg-1];
            nseg --;
        }
    }

    if (DEBUG) {
        std::cout << "segment: ";
        print(std::cout, segment, nseg*NODES_PER_FACET);
        std::cout << '\n';
        std::cout << "segflag: ";
        print(std::cout, segflag, nseg);
        std::cout << '\n';
    }
}


void new_mesh(const Param &param, Variables &var, int bad_quality,
              const array_t &old_coord, const conn_t &old_connectivity,
              const segment_t &old_segment, const segflag_t &old_segflag)
{
    // We don't want to refine large elements during remeshing,
    // so using the domain size as the max area
    double max_elem_size;
#ifdef THREED
    max_elem_size = param.mesh.xlength * param.mesh.ylength * param.mesh.zlength;
#else
    max_elem_size = param.mesh.xlength * param.mesh.zlength;
#endif
    const int vertex_per_polygon = 3;
    const double min_dist = std::pow(param.mesh.smallest_size*sizefactor, 1./NDIMS) * param.mesh.resolution;

    int_vec facet_polygons[6];
    assemble_facet_polygons(var, old_coord, facet_polygons);

    // create a copy of old_coord and old_segment
    double *qcoord = new double[old_coord.num_elements()];
    std::memcpy(qcoord, old_coord.data(), sizeof(double)*old_coord.num_elements());
    int *qsegment = new int[old_segment.num_elements()];
    std::memcpy(qsegment, old_segment.data(), sizeof(int)*old_segment.num_elements());
    int *qsegflag = new int[old_segflag.num_elements()];
    std::memcpy(qsegflag, old_segflag.data(), sizeof(int)*old_segflag.num_elements());

    int old_nnode = old_coord.size();
    int old_nseg = old_segment.size();

    int_vec points_to_delete;

    if (bad_quality == 3) {
        // some elements are too small, delete them
        int_vec tiny_elems;
        find_tiny_element(param, *var.volume, tiny_elems);

        if (tiny_elems.size() > 0) {
            find_points_of_tiny_elem(old_coord, old_connectivity, *var.volume,
                                     tiny_elems, old_nnode, qcoord, *var.bcflag, points_to_delete,
                                     is_boundary);
        }
    }

    // function pointer
    bool (* excl_func)(uint) = NULL;
    switch (param.mesh.remeshing_option) {
    case 0:
        excl_func = &is_boundary;
        break;
    case 1: {
        excl_func = &is_boundary;
        flatten_bottom(*var.bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        break;
    }
    case 2: {
        excl_func = &is_boundary;
        new_bottom(*var.bcflag, qcoord, -param.mesh.zlength,
                   points_to_delete, min_dist, qsegment, qsegflag, old_nseg);
        delete_points(points_to_delete, old_nnode, old_nseg,
                      qcoord, qsegment);
        delete_facets(old_nseg, qsegment, qsegflag);
        break;
    }
    case 10: {
        excl_func = &is_corner;
        // mark points on segment to delete ...
        delete_points_on_boundary(points_to_delete, old_nnode, old_nseg,
                                  qcoord, qsegment, *var.bcflag, min_dist);
        delete_facets(old_nseg, qsegment, qsegflag);
        break;
    }
    case 11: {
        excl_func = &is_corner;
        flatten_bottom(*var.bcflag, qcoord, -param.mesh.zlength,
                       points_to_delete, min_dist);
        delete_points_on_boundary(points_to_delete, old_nnode, old_nseg,
                                  qcoord, qsegment, *var.bcflag, min_dist);
        delete_facets(old_nseg, qsegment, qsegflag);
        break;
    }
    default:
        std::cerr << "Error: unknown remeshing_option: " << param.mesh.remeshing_option << '\n';
        std::exit(1);
    }

    // new mesh
    int new_nnode, new_nelem, new_nseg;
    double *pcoord, *pregattr;
    int *pconnectivity, *psegment, *psegflag;
    Mesh mesh_param = param.mesh;
    points_to_new_mesh(mesh_param, old_nnode, qcoord,
                       old_nseg, qsegment, qsegflag,
                       0, NULL,
                       max_elem_size, vertex_per_polygon,
                       new_nnode, new_nelem, new_nseg,
                       pcoord, pconnectivity, psegment, psegflag, pregattr);

    array_t new_coord(pcoord, new_nnode);
    conn_t new_connectivity(pconnectivity, new_nelem);

    uint_vec new_bcflag(new_nnode);
    create_boundary_flags2(new_bcflag, new_nseg, psegment, psegflag);

    // coarsening
    {
        // deleting (non-boundary) nodes to avoid having tiny elements
        double_vec new_volume(new_nelem);
        compute_volume(new_coord, new_connectivity, new_volume);

        int_vec tiny_elems;
        find_tiny_element(param, new_volume, tiny_elems);

        int_vec points_to_delete;
        if (tiny_elems.size() > 0) {
            find_points_of_tiny_elem(new_coord, new_connectivity, new_volume,
                                     tiny_elems, old_nnode, qcoord, *var.bcflag, points_to_delete,
                                     excl_func);
        }

        // deleting boundary nodes associated with tiny facets
        if (points_to_delete.size() > 0) {
            int q_nnode = old_nnode;
            delete_points_on_boundary(points_to_delete, q_nnode, old_nseg,
                                      qcoord, qsegment, new_bcflag, min_dist);
            delete_facets(old_nseg, qsegment, qsegflag);

            delete [] psegment;
            delete [] psegflag;

            // lessen the quality constraint, so less new points got inserted to the mesh
            mesh_param.min_angle *= 0.9;
            mesh_param.max_ratio *= 1.1;
            points_to_new_mesh(mesh_param, q_nnode, qcoord,
                               old_nseg, qsegment, qsegflag,
                               0, NULL,
                               max_elem_size, vertex_per_polygon,
                               new_nnode, new_nelem, new_nseg,
                               pcoord, pconnectivity, psegment, psegflag, pregattr);

            new_coord.reset(pcoord, new_nnode);
            new_connectivity.reset(pconnectivity, new_nelem);
        }
    }

    var.nnode = new_nnode;
    var.nelem = new_nelem;
    var.nseg = new_nseg;
    var.coord->steal_ref(new_coord);
    var.connectivity->steal_ref(new_connectivity);
    var.segment->reset(psegment, var.nseg);
    var.segflag->reset(psegflag, var.nseg);

    delete [] qcoord;
    delete [] qsegment;
    delete [] qsegflag;
}

} // anonymous namespace


int bad_mesh_quality(const Param &param, const Variables &var, int &index)
{
    /* Check the quality of the mesh, return 0 if the mesh quality (by several
     * measures) is good. Non-zero returned values indicate --
     * 1: an element has bad quality (too acute / narrow / flat).
     * 2: a bottom node has moved too far away from the flat bottom.
     * 3: an element is smaller than (mesh.smallest_size * [volume of a equilateral triangle/tetrahedron
     *    of side = mesh.resolution]).
     */

    // check tiny elements
    const double smallest_vol = param.mesh.smallest_size * sizefactor * std::pow(param.mesh.resolution, NDIMS);
    for (int e=0; e<var.nelem; e++) {
        if ((*var.volume)[e] < smallest_vol) {
            index = e;
            std::cout << "    The size of element #" << index << " is too small.\n";
            return 3;
        }
    }

    // check if any bottom node is too far away from the bottom depth
    if (param.mesh.remeshing_option == 1 ||
        param.mesh.remeshing_option == 2 ||
        param.mesh.remeshing_option == 11) {
        double bottom = - param.mesh.zlength;
        const double dist_ratio = 0.25;
        for (int i=0; i<var.nnode; ++i) {
            if (is_bottom((*var.bcflag)[i])) {
                double z = (*var.coord)[i][NDIMS-1];
                if (std::fabs(z - bottom) > dist_ratio * param.mesh.resolution) {
                    index = i;
                    std::cout << "    Node #" << i << " is too far from the bottm: z = " << z << "\n";
                    return 2;
                }
            }
        }
    }

    // check element distortion
    int worst_elem;
    double q = worst_elem_quality(*var.coord, *var.connectivity,
                                  *var.volume, *var.elquality, worst_elem);
#ifdef THREED
    // normalizing q so that its magnitude is about the same in 2D and 3D
    q = std::pow(q, 1.0/3);
#endif
    if (q < param.mesh.min_quality) {
        index = worst_elem;
        std::cout << "    Element #" << worst_elem << " has mesh quality = " << q << ".\n";
        return 1;
    }
    return 0;
}


void remesh(const Param &param, Variables &var, int bad_quality)
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

        new_mesh(param, var, bad_quality, old_coord, old_connectivity,
                 old_segment, old_segflag);

        // interpolating fields defined on elements
        nearest_neighbor_interpolation(var, old_coord, old_connectivity);

        // interpolating fields defined on nodes
        barycentric_node_interpolation(var, old_coord, old_connectivity);

        // remap markers. elemmarkers are updated here, too.
        remap_markers(param, var, old_coord, old_connectivity);
  
        // old_coord et al. are destroyed before exiting this block
    }

    // memory for new fields
    reallocate_variables(param, var);

    // updating other arrays
    delete var.bcflag;
    create_boundary_flags(var);
    for (int i=0; i<6; ++i) {
        var.bnodes[i].clear();
        var.bfacets[i].clear();
    }
    create_boundary_nodes(var);
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

    if (param.sim.output_during_remeshing) {
        // the following variables need to be re-computed only when we are
        // outputing right after remeshing
        update_strain_rate(var, *var.strain_rate);
        update_force(param, var, *var.force);
        int junk;
        worst_elem_quality(*var.coord, *var.connectivity,
                           *var.volume, *var.elquality, junk);
    }

    std::cout << "  Remeshing finished.\n";
}


