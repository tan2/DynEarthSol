#include "iostream"

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
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


static void barycentric_node_interpolation(Variables &var, const array_t &old_coord,
                                           const conn_t &old_connectivity)
{
    // for each new coord point, find the enclosing old element
    Barycentric_transformation bary(old_coord, old_connectivity);
    int_vec el(var.nnode);
    Array2D<double,NODES_PER_ELEM> brc_coord(var.nnode);
    {
        // ANN requires double** as input
        double **points = new double*[old_coord.size()];
        for (int i=0; i<var.nnode; i++) {
            points[i] = const_cast<double*>(old_coord[i]);
        }
        ANNkd_tree kdtree(points, old_coord.size(), NDIMS);

        const std::vector<int_vec> &old_support = *var.support;

        const int k = 1;
        const double eps = 0;
        int *nn_idx = new int[k];
        double *dd = new double[k];

        for (int i=0; i<var.nnode; i++) {
            double *q = (*var.coord)[i];
            // find the nearest point nn in old_coord
            kdtree.annkSearch(q, k, nn_idx, dd, eps);
            int nn = nn_idx[0];

            // loop over (old) elements surrounding nn to find
            // the element that is enclosing q
            int e;
            double r[NDIMS];
            const int_vec &nn_elem = old_support[nn];
            for (int j=0; j<nn_elem.size(); j++) {
                e = nn_elem[j];
                bary.transform(q, e, r);
                if (bary.is_inside(r)) {
                    goto found;
                }
            }
        not_found:
            // q is outside the old domain
            // using nearest old_coord instead
            bary.transform(points[nn], e, r);

        found:
            el[i] = e;
            double sum = 0;
            for (int d=0; d<NDIMS; d++) {
                brc_coord[i][d] = r[d];
                sum += r[d];
            }
            brc_coord[i][NODES_PER_ELEM-1] = 1 - sum;
        }

        delete [] nn_idx;
        delete [] dd;
        // XXX: this line causes double-free error on run-time, why?
        //delete [] points;

        // print(std::cout, *var.coord);
        // std::cout << '\n';
        // print(std::cout, el);
        // std::cout << '\n';
        // print(std::cout, bar);
        // std::cout << '\n';
    }

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

    // memory for new fields
    reallocate_variables(param, var);

    // interpolating fields
    nearest_neighbor_interpolation(var, old_coord, old_connectivity);
    barycentric_node_interpolation(var, old_coord, old_connectivity);

    // arrays of new mesh

    // updating other arrays

    free(old_coord.data());
    free(old_connectivity.data());
    free(old_segment.data());
    free(old_segflag.data());

    std::cout << "  Remeshing finished.\n";
}


