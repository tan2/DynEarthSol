#include "iostream"

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "mesh.hpp"
#include "nn-interpolation.hpp"


void find_nearest_neighbor(Variables &var, const array_t &old_coord,
                           const conn_t &old_connectivity, int_vec &idx)
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
    // Note: kdtree.annkSearch() is not thread-safe, cannot use openmp in this loop
    for(int e=0; e<var.nelem; e++) {
        double *q = new_center[e];
        kdtree.annkSearch(q, k, nn_idx, dd, eps);
        idx[e] = nn_idx[0];
    }

    delete [] nn_idx;
    delete [] dd;
    delete [] new_center[0];
    delete [] new_center;
    delete [] old_center[0];
    delete [] old_center;
}


static void inject_field(const int_vec &idx, const double_vec &source, double_vec &target)
{
    #pragma omp parallel for default(none)          \
        shared(idx, source, target)
    for (std::size_t i=0; i<target.size(); i++) {
        int n = idx[i];
        target[i] = source[n];
    }
}


static void inject_field(const int_vec &idx, const tensor_t &source, tensor_t &target)
{
    #pragma omp parallel for default(none)          \
        shared(idx, source, target)
    for (std::size_t i=0; i<target.size(); i++) {
        int n = idx[i];
        for (int d=0; d<NSTR; d++) {
            target[i][d] = source[n][d];
        }
    }
}


static void nn_interpolate_elem_fields(Variables &var, const int_vec &idx)
{
    const int n = var.nnode;
    const int e = var.nelem;

    double_vec *a = new double_vec(e);

    inject_field(idx, *var.plstrain, *a);
    delete var.plstrain;
    var.plstrain = a;

    tensor_t *b;
    b = new tensor_t(e);
    inject_field(idx, *var.strain, *b);
    delete var.strain;
    var.strain = b;

    b = new tensor_t(e);
    inject_field(idx, *var.stress, *b);
    delete var.stress;
    var.stress = b;

}


void nearest_neighbor_interpolation(Variables &var, const array_t &old_coord,
                                    const conn_t &old_connectivity)
{
    int_vec idx(var.nelem);
    find_nearest_neighbor(var, old_coord, old_connectivity, idx);

    nn_interpolate_elem_fields(var, idx);

}
