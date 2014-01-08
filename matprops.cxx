
#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "utils.hpp"
#include "matprops.hpp"


namespace {

    double arithmetic_mean(const double_vec &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        for (std::size_t i=0; i<s.size(); i++) {
            result += n[i] * s[i];
            m += n[i];
        }
        return result / m;
    }


    double harmonic_mean(const double_vec &s, const int_vec &n)
    {
        if (s.size() == 1) return s[0];

        double result = 0;
        int m = 0;
        for (std::size_t i=0; i<s.size(); i++) {
            result += n[i] / s[i];
            m += n[i];
        }
        return m / result;
    }

}


MatProps::MatProps(const Param& p, const Variables& var) :
  rheol_type(p.mat.rheol_type),
  nmat(p.mat.nmat),
  visc_min(p.mat.visc_min),
  visc_max(p.mat.visc_max),
  tension_max(p.mat.tension_max),
  therm_diff_max(p.mat.therm_diff_max),
  rho0(p.mat.rho0),
  alpha(p.mat.alpha),
  bulk_modulus(p.mat.bulk_modulus),
  shear_modulus(p.mat.shear_modulus),
  visc_exponent(p.mat.visc_exponent),
  visc_coefficient(p.mat.visc_coefficient),
  visc_activation_energy(p.mat.visc_activation_energy),
  heat_capacity(p.mat.heat_capacity),
  therm_cond(p.mat.therm_cond),
  pls0(p.mat.pls0),
  pls1(p.mat.pls1),
  cohesion0(p.mat.cohesion0),
  cohesion1(p.mat.cohesion1),
  friction_angle0(p.mat.friction_angle0),
  friction_angle1(p.mat.friction_angle1),
  dilation_angle0(p.mat.dilation_angle0),
  dilation_angle1(p.mat.dilation_angle1),
  coord(*var.coord),
  connectivity(*var.connectivity),
  temperature(*var.temperature),
  stress(*var.stress),
  strain_rate(*var.strain_rate),
  elemmarkers(*var.elemmarkers)
{}


MatProps::~MatProps()
{}


double MatProps::bulkm(int e) const
{
    return harmonic_mean(bulk_modulus, elemmarkers[e]);
}


double MatProps::shearm(int e) const
{
    return harmonic_mean(shear_modulus, elemmarkers[e]);
}


double MatProps::visc(int e) const
{
    const double gas_constant = 8.3144;
    const double min_strain_rate = 1e-30;

    // average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    // strain-rate
    double edot = second_invariant(strain_rate[e]);
    // min strain rate to prevent viscosity -> inf
    edot = std::max(edot, min_strain_rate);

    // viscosity law from Chen and Morgan, JGR, 1990
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        double pow = 1 / visc_exponent[m] - 1;
        double pow1 = -1 / visc_exponent[m];
        double visc0 = 0.25 * std::pow(edot, pow) * std::pow(0.75 * visc_coefficient[m], pow1)
            * std::exp(visc_activation_energy[m] / (visc_exponent[m] * gas_constant * T)) * 1e6;
        result += elemmarkers[e][m] * visc0;
        n += elemmarkers[e][m];
    }

    double visc = result / n;

    // applying min & max limits
    visc = std::min(std::max(visc, visc_min), visc_max);

    return visc;
}


void MatProps::plastic_weakening(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening) const
{
    // TODO: compute average plastic properties
    const int mat = 0;
    double c, f, d, h;
    if (pls <= pls0[mat]) {
        // no weakening yet
        c = cohesion0[mat];
        f = friction_angle0[mat];
        d = dilation_angle0[mat];
        h = 0;
    }
    else if (pls < pls1[mat]) {
        // linear weakening
        double p = (pls - pls0[mat]) / (pls1[mat] - pls0[mat]);
        c = cohesion0[mat] + p * (cohesion1[mat] - cohesion0[mat]);
        f = friction_angle0[mat] + p * (friction_angle1[mat] - friction_angle0[mat]);
        d = dilation_angle0[mat] + p * (dilation_angle1[mat] - dilation_angle0[mat]);
        h = (cohesion1[mat] - cohesion0[mat]) / (pls1[mat] - pls0[mat]);
    }
    else {
        // saturated weakening
        c = cohesion1[mat];
        f = friction_angle1[mat];
        d = dilation_angle1[mat];
        h = 0;
    }

    cohesion = c;
    friction_angle = f;
    dilation_angle = d;
    hardening = h;
}


void MatProps::plastic_props(int e, double pls,
                             double& amc, double& anphi, double& anpsi,
                             double& hardn, double& ten_max) const
{
    double cohesion, phi, psi;

    plastic_weakening(e, pls, cohesion, phi, psi, hardn);

    // derived variables
    double sphi = std::sin(phi * DEG2RAD);
    double spsi = std::sin(psi * DEG2RAD);

    anphi = (1 + sphi) / (1 - sphi);
    anpsi = (1 + spsi) / (1 - spsi);
    amc = 2 * cohesion * std::sqrt(anphi);

    ten_max = (phi == 0)? tension_max : std::min(tension_max, cohesion/std::tan(phi*DEG2RAD));
}


double MatProps::rho(int e) const
{
    const double celsius0 = 273;

    // average temperature of this element
    double T = 0;
    const int *conn = connectivity[e];
    for (int i=0; i<NODES_PER_ELEM; ++i) {
        T += temperature[conn[i]];
    }
    T /= NODES_PER_ELEM;

    double TinCelsius = T - celsius0;
    double result = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        // TODO: compressibility
        result += (rho0[m] - alpha[m] * TinCelsius) * elemmarkers[e][m];
        n += elemmarkers[e][m];
    }
    return result / n;
}


double MatProps::cp(int e) const
{
    return arithmetic_mean(heat_capacity, elemmarkers[e]);
}


double MatProps::k(int e) const
{
    return arithmetic_mean(therm_cond, elemmarkers[e]);
}


