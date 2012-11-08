#ifndef DYNEARTHSOL3D_MATPROPS_HPP
#define DYNEARTHSOL3D_MATPROPS_HPP

#include "parameters.hpp"

class MatProps
{
public:
    MatProps(int nmat, int rheol_type);

    const int nmat;
    const int rheol_type;

    double bulkm(int e) const;
    double shearm(int e) const;
    double visc(int e) const;

    double density(int e) const;
    double cp(int e) const;
    double k(int e) const;

    void plastic_props(int e, double pls,
                       double& amc, double& anphi, double& anpsi,
                       double& hardn, double& ten_max) const;

    // XXX: take values from param
    const double visc_min = 1e21;
    const double visc_max = 1e21;
    const double therm_diff_max = 1e-6;

    const static int rh_elastic = 1 << 0;
    const static int rh_viscous = 1 << 1;
    const static int rh_plastic = 1 << 2;
    const static int rh_maxwell = rh_elastic | rh_viscous;
    const static int rh_ep = rh_elastic | rh_plastic;
    const static int rh_evp = rh_elastic | rh_viscous | rh_plastic;

private:
    double_vec bulk_modulus, shear_modulus;
    double_vec rho0, alpha;
    double_vec heat_capacity, therm_cond;
    double ten_off;
};


#endif
