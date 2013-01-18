#include "algorithm"
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

            // std::cout << i << " ";
            // print(std::cout, q, NDIMS);
            // std::cout << " " << nn << " " << dd[0] << '\n';

            // loop over (old) elements surrounding nn to find
            // the element that is enclosing q
            int e;
            double r[NDIMS];
            const int_vec &nn_elem = old_support[nn];
            for (int j=0; j<nn_elem.size(); j++) {
                e = nn_elem[j];
                bary.transform(q, e, r);
                if (bary.is_inside(r)) {
                    // std::cout << e << " ";
                    // print(std::cout, r, NDIMS);
                    // std::cout << '\n';
                    goto found;
                }
            }
        not_found:
            {
                // Situation: q is in the upper element, but its nearest point is o!
                // we won't find the enclosing element with the method above
                //     x
                //    / \   <-- this is a large triangle
                //   / q \
                //  x---- x
                //   \-o-/   <-- this is a small triangle
                //

                // this array contains the elements that have been searched so far
                int_vec searched;
                for (int j=0; j<nn_elem.size(); j++) {
                    searched.push_back(nn_elem[j]);
                }
                // print(std::cout, searched);
                // std::cout << " ... \n";

                // search through elements that are neighbors of nn_elem
                for (int j=0; j<nn_elem.size(); j++) {
                    int ee = nn_elem[j];
                    const int *conn = old_connectivity[ee];
                    for (int m=0; m<NODES_PER_ELEM; m++) {
                        // np is a node close to q
                        int np = conn[m];
                        const int_vec &np_elem = old_support[np];
                        for (int j=0; j<np_elem.size(); j++) {
                            e = np_elem[j];
                            auto it = std::find(searched.begin(), searched.end(), e);
                            if (it != searched.end()) {
                                // this element has been searched before
                                continue;
                            }
                            searched.push_back(e);
                            bary.transform(q, e, r);
                            // std::cout << e << " ";
                            // print(std::cout, r, NDIMS);
                            // std::cout << " ... \n";
                            if (bary.is_inside(r)) {
                                goto found;
                            }
                        }
                    }
                }
            }
            {
                std::cout << "New node is outside of the old domain. \n";

                // Situation: q must be outside the old domain
                // using nearest old_coord instead
                bary.transform(points[nn], nn_elem[0], r);
            }
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


