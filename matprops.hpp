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
    double density(int e) const;
    double visc(int e) const;

    const static int rh_elastic = 0x00000001;
    const static int rh_viscous = 0x00000010;
    const static int rh_plastic = 0x00000100;
    const static int rh_maxwell = rh_elastic & rh_viscous;
    const static int rh_ep = rh_elastic & rh_plastic;
    const static int rh_evp = rh_elastic & rh_viscous & rh_plastic;

private:
    double_vec bulk_modulus, shear_modulus;
    double_vec rho0, alpha;
    double_vec heat_capacity, therm_cond;
};


#endif
