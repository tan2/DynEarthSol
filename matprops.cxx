
#include <algorithm>

#include "parameters.hpp"
#include "matprops.hpp"

void create_matprops(const Param &par, Variables &var)
{
    var.mat = new MatProps(1, MatProps::rh_evp);
}


MatProps::MatProps(int n, int rh) :
    nmat(n), rheol_type(rh)
{
    bulk_modulus.reserve(n);
    shear_modulus.reserve(n);

    rho0.reserve(n);
    alpha.reserve(n);
    heat_capacity.reserve(n);
    therm_cond.reserve(n);

    std::fill_n(bulk_modulus.begin(), n, 128.2e9);
    std::fill_n(shear_modulus.begin(), n, 80.5e9);

    std::fill_n(rho0.begin(), n, 3210);
    std::fill_n(alpha.begin(), n, 3e-5);
    std::fill_n(heat_capacity.begin(), n,1000);
    std::fill_n(therm_cond.begin(), n, 3);
}
