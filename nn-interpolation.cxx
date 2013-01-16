#include "iostream"

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"

#include "nn-interpolation.hpp"


double** elem_center(const array_t &coord, const conn_t &connectivity)
{
    /* Returns the centroid of the elements.
     * Note: center[0] == tmp
     * The caller is responsible to delete [] center[0] and center!
     */
    int nelem = connectivity.size();
    double *tmp = new double[nelem*NDIMS];
    double **center = new double*[nelem];
    for(int e=0; e<nelem; e++) {
        const int* conn = connectivity[e];
        center[e] = tmp + e*NDIMS;
        for(int d=0; d<NDIMS; d++) {
            double sum = 0;
            for(int k=0; k<NODES_PER_ELEM; k++) {
                sum += coord[conn[k]][d];
            }
            center[e][d] = sum / NODES_PER_ELEM;
        }
    }
    return center;
}


void nearest_neighbor_interpolation(Variables &var, const array_t &old_coord,
                                    const conn_t &old_connectivity)
{
    std::cout << "Constructing a kd-tree.\n";
    // kdtree requires the coordinate as double**
    double **old_center = elem_center(old_coord, old_connectivity);
    ANNkd_tree kdtree(old_center, old_connectivity.size(), NDIMS);

    std::cout << "Searching nearest neighbor in the kd-tree.\n";
    double **new_center = elem_center(*var.coord, *var.connectivity);
    const int k = 1;
    const double eps = 0;
    int *nn_idx = new int[k];
    double *dd = new double[k];
    for(int i=0; i<var.nnode; i++) {
        double *q = new_center[i];
        kdtree.annkSearch(q, k, nn_idx, dd, eps);
    }

    delete [] nn_idx;
    delete [] dd;
    delete [] new_center[0];
    delete [] new_center;
    delete [] old_center[0];
    delete [] old_center;
}
