#include "iostream"

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"

#include "nn-interpolation.hpp"


void nearest_neighbor_interpolation(Variables &var, const array_t &old_coord)
{
    std::cout << "Constructing a kd-tree.\n";
    // kdtree requires the coordinate as double**
    double **pa = new double*[var.nnode];
    for(int i=0; i<var.nnode; i++) {
        pa[i] = (*var.coord)[i];
    }
    ANNkd_tree kdtree(pa, old_coord.size(), NDIMS);

    std::cout << "Searching nearest neighbor in the kd-tree.\n";
    const int k = 1;
    const double eps = 0;
    int *nn_idx = new int[k];
    double *dd = new double[k];
    for(int i=0; i<var.nnode; i++) {
        double *q = (*var.coord)[i];
        kdtree.annkSearch(q, k, nn_idx, dd, eps);
    }
}
