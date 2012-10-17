
#include <algorithm>

#include "matprops.hpp"

MatProps::MatProps(int n, int rh) :
    nmat(n), rheol_type(rh)
{
    bulk_modulus.reserve(n);
    shear_modulus.reserve(n);

    rho0.reserve(n);
    alpha.reserve(n);
    heat_capacity.reserve(n);
    therm_cond.reserve(n);

    // TODO: get material properties from cfg file
    std::fill_n(bulk_modulus.begin(), n, 128.2e9);
    std::fill_n(shear_modulus.begin(), n, 80.5e9);

    std::fill_n(rho0.begin(), n, 3210);
    std::fill_n(alpha.begin(), n, 3e-5);
    std::fill_n(heat_capacity.begin(), n,1000);
    std::fill_n(therm_cond.begin(), n, 3);
}


double MatProps::bulkm(int e) const
{
    // TODO: compute average bulk modulus
    return bulk_modulus[0];
}


double MatProps::shearm(int e) const
{
    // TODO: compute average shear modulus
    return shear_modulus[0];
}


double MatProps::density(int e) const
{
    // TODO: compute average density with thermal expansion
    return rho0[0];
}


double MatProps::visc(int e) const
{
    // TODO: compute average viscosity
    return 1e20;
}
