
#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "matprops.hpp"


MatProps::MatProps(const Param& p) :
  rheol_type(p.mat.rheol_type),
  nmat(p.mat.nmat),
  visc_min(p.mat.visc_min),
  visc_max(p.mat.visc_max),
  therm_diff_max(p.mat.therm_diff_max),
  rho0(p.mat.rho0),
  alpha(p.mat.alpha),
  bulk_modulus(p.mat.bulk_modulus),
  shear_modulus(p.mat.shear_modulus),
  heat_capacity(p.mat.heat_capacity),
  therm_cond(p.mat.therm_cond)
{
    // TODO: get material properties from cfg file
    ten_off = 1e9;
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


double MatProps::visc(int e) const
{
    // TODO: compute average viscosity
    return visc_max;
}


void MatProps::plastic_props(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max) const
{
    // TODO: compute average plastic properties

    // plastic properties due to strain weakening
    double pls_seg[2] = {0.0, 0.1};
    double coh_seg[2] = {4e7, 4e5};  // in Pa
    double fric_seg[2] = {15, 1};  // in degree
    double dilat_seg[2] = {0, 0};  // in degree

    double c, f, d, h;

    if (pls < pls_seg[0]) {
        // no weakening yet
        c = coh_seg[0];
        f = fric_seg[0];
        d = dilat_seg[0];
        h = 0;
    }
    else if (pls < pls_seg[1]) {
        // linear weakening
        double p = (pls - pls_seg[0]) / (pls_seg[1] - pls_seg[0]);
        c = coh_seg[0] + p * (coh_seg[1] - coh_seg[0]);
        f = fric_seg[0] + p * (fric_seg[1] - fric_seg[0]);
        d = dilat_seg[0] + p * (dilat_seg[1] - dilat_seg[0]);
        h = (coh_seg[1] - coh_seg[0]) / (pls_seg[1] - pls_seg[0]);
    }
    else {
        // saturated weakening
        c = coh_seg[1];
        f = fric_seg[1];
        d = dilat_seg[1];
        h = 0;
    }

    hardn = h;

    double cohesion = c;
    double phi = f;
    double psi = d;

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? ten_off : std::min(ten_off, cohesion/std::tan(phi*DEG2RAD));
}


double MatProps::density(int e) const
{
    // TODO: compute average density with thermal expansion
    return rho0[0];
}


double MatProps::cp(int e) const
{
    // TODO: compute average heat capacity
    return heat_capacity[0];
}


double MatProps::k(int e) const
{
    // TODO: compute average thermal conductivity
    return therm_cond[0];
}


