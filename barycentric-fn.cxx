
#include "barycentric-fn.hpp"


Barycentric_transformation::Barycentric_transformation(const array_t &coord, const conn_t &connectivity)
    : N(connectivity.size()),
      coeff_(connectivity.size())
{
    #pragma omp parallel for default(none) \
        shared(coord, connectivity)
    for (int e=0; e<N; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

#ifdef THREED
        int n3 = connectivity[e][3];
        const double *d = coord[n3];

        compute_coeff3d(a, b, c, d, coeff_[e]);
#else
        compute_coeff2d(a, b, c, coeff_[e]);
#endif
    }
}


Barycentric_transformation::~Barycentric_transformation() {};


void Barycentric_transformation::transform(const double *point, int e, double *result) const
{
    const double *cf = coeff_[e];
    for (int d=0; d<NDIMS; d++) {
        result[d] = cf[index(0,d)];
        for (int i=0; i<NDIMS; i++)
            result[d] += cf[index(i+1,d)]*point[i];
    }
}


bool Barycentric_transformation::is_inside_elem(const double *point, int elem) const
{
    double r[NDIMS];
    transform(point, elem, r);
    return is_inside(r);
}


bool Barycentric_transformation::is_inside(const double *r) const
{
    const double tolerance = 1e-14;

#ifdef THREED
    if (r[0] >= -tolerance &&
        r[1] >= -tolerance &&
        r[2] >= -tolerance &&
        (r[0] + r[1] + r[2]) <= 1 + tolerance)
        return 1;
#else
    if (r[0] >= -tolerance &&
        r[1] >= -tolerance &&
        (r[0] + r[1]) <= 1 + tolerance)
        return 1;
#endif
    return 0;
}


inline int Barycentric_transformation::index(int node, int dim) const
{
    return node*NDIMS + dim;
}


#ifdef THREED

void Barycentric_transformation::compute_coeff3d(const double *a,
                                                 const double *b,
                                                 const double *c,
                                                 const double *d,
                                                 double *coeff_e)
{
    //TODO
}

#else

void Barycentric_transformation::compute_coeff2d(const double *a,
                                                 const double *b,
                                                 const double *c,
                                                 double *coeff_e)
{
    double det = (a[0]*b[1] - a[1]*b[0] +
                  b[0]*c[1] - b[1]*c[0] +
                  c[0]*a[1] - c[1]*a[0]);

    coeff_e[index(0,0)] = (b[0]*c[1] - b[1]*c[0]) / det;
    coeff_e[index(0,1)] = (c[0]*a[1] - c[1]*a[0]) / det;
    coeff_e[index(1,0)] = (b[1] - c[1]) / det;
    coeff_e[index(1,1)] = (c[1] - a[1]) / det;
    coeff_e[index(2,0)] = (c[0] - b[0]) / det;
    coeff_e[index(2,1)] = (a[0] - c[0]) / det;
}

#endif


