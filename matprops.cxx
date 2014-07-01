
#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "utils.hpp"
#include "matprops.hpp"


namespace {

    double arithmetic_mean(const VectorBase &s, const int_vec &n)
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


    double harmonic_mean(const VectorBase &s, const int_vec &n)
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


    class Vector : public VectorBase {
    private:
        const double_vec &a;
        const int len;
    public:
        Vector(const double_vec &a_, int len_) :
            a(a_), len(len_)
        {}

        double operator[](std::size_t i) const
        {
            return a[i];
        }

        std::size_t size() const
        {
            return a.size();
        }
    };


    class Vector1 : public VectorBase {
    private:
        const double d;
        const int len;
    public:
        Vector1(const double_vec &a_, int len_) :
            d(a_[0]), len(len)
        {}

        double operator[](std::size_t i) const
        {
            return d; // always return the same element
        }

        std::size_t size() const
        {
            return 1;
        }
    };


}


VectorBase* VectorBase::create(const double_vec &a, int len)
{
    if (a.size() == len)
        return new Vector(a, len);
    if (a.size() == 1 && len != 1)
        return new Vector1(a, len);

    std::cerr << "Error: incorrect parameters received in VectorBase::create() at "
              << __FILE__ << ':' << __LINE__ << '\n';
    std::exit(12);
    return NULL;
}


MatProps::MatProps(const Param& p, const Variables& var) :
  rheol_type(p.mat.rheol_type),
  nmat(p.mat.nmat),
  visc_min(p.mat.visc_min),
  visc_max(p.mat.visc_max),
  tension_max(p.mat.tension_max),
  therm_diff_max(p.mat.therm_diff_max),
  coord(*var.coord),
  connectivity(*var.connectivity),
  temperature(*var.temperature),
  stress(*var.stress),
  strain_rate(*var.strain_rate),
  elemmarkers(*var.elemmarkers)
{
    rho0 = VectorBase::create(p.mat.rho0, nmat);
    alpha = VectorBase::create(p.mat.alpha, nmat);
    bulk_modulus = VectorBase::create(p.mat.bulk_modulus, nmat);
    shear_modulus = VectorBase::create(p.mat.shear_modulus, nmat);
    visc_exponent = VectorBase::create(p.mat.visc_exponent, nmat);
    visc_coefficient = VectorBase::create(p.mat.visc_coefficient, nmat);
    visc_activation_energy = VectorBase::create(p.mat.visc_activation_energy, nmat);
    heat_capacity = VectorBase::create(p.mat.heat_capacity, nmat);
    therm_cond = VectorBase::create(p.mat.therm_cond, nmat);
    pls0 = VectorBase::create(p.mat.pls0, nmat);
    pls1 = VectorBase::create(p.mat.pls1, nmat);
    cohesion0 = VectorBase::create(p.mat.cohesion0, nmat);
    cohesion1 = VectorBase::create(p.mat.cohesion1, nmat);
    friction_angle0 = VectorBase::create(p.mat.friction_angle0, nmat);
    friction_angle1 = VectorBase::create(p.mat.friction_angle1, nmat);
    dilation_angle0 = VectorBase::create(p.mat.dilation_angle0, nmat);
    dilation_angle1 = VectorBase::create(p.mat.dilation_angle1, nmat);
}


MatProps::~MatProps()
{
    delete rho0;
    delete alpha;
    delete bulk_modulus;
    delete shear_modulus;
    delete visc_exponent;
    delete visc_coefficient;
    delete visc_activation_energy;
    delete heat_capacity;
    delete therm_cond;
    delete pls0;
    delete pls1;
    delete cohesion0;
    delete cohesion1;
    delete friction_angle0;
    delete friction_angle1;
    delete dilation_angle0;
    delete dilation_angle1;
}


double MatProps::bulkm(int e) const
{
    return harmonic_mean(*bulk_modulus, elemmarkers[e]);
}


double MatProps::shearm(int e) const
{
    return harmonic_mean(*shear_modulus, elemmarkers[e]);
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
        double pow = 1 / (*visc_exponent)[m] - 1;
        double pow1 = -1 / (*visc_exponent)[m];
        double visc0 = 0.25 * std::pow(edot, pow) * std::pow(0.75 * (*visc_coefficient)[m], pow1)
            * std::exp((*visc_activation_energy)[m] / ((*visc_exponent)[m] * gas_constant * T)) * 1e6;
        result += elemmarkers[e][m] / visc0;
        n += elemmarkers[e][m];
    }

    double visc = n / result;

    // applying min & max limits
    visc = std::min(std::max(visc, visc_min), visc_max);

    return visc;
}


void MatProps::plastic_weakening(int e, double pls,
                                 double &cohesion, double &friction_angle,
                                 double &dilation_angle, double &hardening) const
{
    double c, f, d, h;
    c = f = d = h = 0;
    int n = 0;
    for (int m=0; m<nmat; m++) {
        int k = elemmarkers[e][m];
        if (k == 0) continue;
        n += k;
        if (pls < (*pls0)[m]) {
            // no weakening yet
            c += (*cohesion0)[m] * k;
            f += (*friction_angle0)[m] * k;
            d += (*dilation_angle0)[m] * k;
            h += 0;
        }
        else if (pls < (*pls1)[m]) {
            // linear weakening
            double p = (pls - (*pls0)[m]) / ((*pls1)[m] - (*pls0)[m]);
            c += ((*cohesion0)[m] + p * ((*cohesion1)[m] - (*cohesion0)[m])) * k;
            f += ((*friction_angle0)[m] + p * ((*friction_angle1)[m] - (*friction_angle0)[m])) * k;
            d += ((*dilation_angle0)[m] + p * ((*dilation_angle1)[m] -(* dilation_angle0)[m])) * k;
            h += ((*cohesion1)[m] - (*cohesion0)[m]) / ((*pls1)[m] - (*pls0)[m]) * k;
        }
        else {
            // saturated weakening
            c += (*cohesion1)[m] * k;
            f += (*friction_angle1)[m] * k;
            d += (*dilation_angle1)[m] * k;
            h += 0;
        }
    }
    cohesion = c / n;
    friction_angle = f / n;
    dilation_angle = d / n;
    hardening = h / n;
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
        result += ((*rho0)[m] - (*alpha)[m] * TinCelsius) * elemmarkers[e][m];
        n += elemmarkers[e][m];
    }
    return result / n;
}


double MatProps::cp(int e) const
{
    return arithmetic_mean(*heat_capacity, elemmarkers[e]);
}


double MatProps::k(int e) const
{
    return arithmetic_mean(*therm_cond, elemmarkers[e]);
}


