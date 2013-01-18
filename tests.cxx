#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
#include "utils.hpp"


void test_barycentric_transformation(Variables &var)
{
    Barycentric_transformation bary(*var.coord, *var.connectivity);

    double p[NDIMS];

    // p is the mid point of element 1
    int e = 1;
    const int *conn = (*var.connectivity)[e];
    for (int d=0; d<NDIMS; d++)
        p[d] = (*var.coord)[conn[0]][d] / NODES_PER_ELEM;
    for (int i=1; i<NODES_PER_ELEM; i++)
        for (int d=0; d<NDIMS; d++)
            p[d] += (*var.coord)[conn[i]][d] / NODES_PER_ELEM;

    double q[NDIMS];
    bary.transform(p, e, q);
    print(std::cout, p, NDIMS);
    // should print out "[0.333333, 0.333333]"
    print(std::cout, q, NDIMS);
}


