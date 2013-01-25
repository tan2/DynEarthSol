#include <algorithm>

#include "constants.hpp"
#include "parameters.hpp"

#include "array2d.hpp"
#include "barycentric-fn.hpp"
#include "sortindex.hpp"
#include "utils.hpp"


void test_array2d()
{

    Array2D<int,3> a(10);
    {
        int *r = new int[10*3];
        Array2D<int,3> a(r, 10);
    }
    for (int i=0; i<10; i++) {
        int *b = a[i];
        for (int j=0; j<3; j++)
            b[j] = i+j;
    }

    int *c = a.data();
    for (int i=0; i<10*3; i++) {
        std::cout << *(c+i) << '\n';
    }
}


void test_barycentric_transformation(Variables &var)
{
    Barycentric_transformation bary(*var.coord, *var.connectivity, *var.volume);

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
    // should print out "[0.333333, 0.333333]" in 2D
    // or "[0.25, 0.25, 0.25]" in 3D
    print(std::cout, q, NDIMS);
    std::cout << '\n';
}


void test_sortindex()
{
    std::vector<std::size_t> list(10), idx(10);
    std::iota(list.begin(), list.begin()+5, 10);
    std::iota(list.begin()+5, list.end(), 0);

    sortindex(list, idx);

    std::cout << "result ...\n";
    for(int i=0; i<10; ++i) {
        std::cout << list[i] << ", " << idx[i] << '\n';
    }

    std::cout << "sorted ...\n";
    for(int i=0; i<10; ++i) {
        std::cout << i << ": " << list[idx[i]] << '\n';
    }
}
