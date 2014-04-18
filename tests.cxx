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


#include "3x3-C/dsyevh3.h"
void test_eigen_decomposition()
{
    // s is a 3x3 tensor, only the upper part is needed.
    double a[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    a[0][0] = 3e4;
    a[1][1] = -1;
    a[2][2] = 3;
    a[0][1] = 2;
    a[0][2] = 4;
    a[1][2] = 2;
    // a[1][0] = a[0][1];
    // a[2][0] = a[0][2];
    // a[2][1] = a[1][2];

    std::cout << "\noriginal matrix:";
    print(std::cout, a[0], 3);
    print(std::cout, a[1], 3);
    print(std::cout, a[2], 3);
    std::cout << "\n";

    // p is eigenvalue, v is eigenvector
    double p[3], v[3][3];

    dsyevh3(a, v, p);

    std::cout << "eigenvalue:";
    print(std::cout, p, 3);
    std::cout << "\neigenvector:";
    print(std::cout, v[0], 3);
    print(std::cout, v[1], 3);
    print(std::cout, v[2], 3);
    std::cout << "\n";

    double ss[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int m=0; m<3; m++) {
        //for(int n=0; n<3; n++) { // using this loop to recover whole matrix
        for(int n=m; n<3; n++) { // using this loop to recover upper half matrix
            for(int k=0; k<3; k++) {
                ss[m][n] += v[m][k] * v[n][k] * p[k];
            }
        }
    }
    std::cout << "\nmatrix:";
    print(std::cout, ss[0], 3);
    print(std::cout, ss[1], 3);
    print(std::cout, ss[2], 3);
    std::cout << "\n";

}


