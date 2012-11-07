#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static void elastic(double bulkm, double shearm, const double* de, double* s)
{
    /* increment the stress s according to the incremental strain de */
    double a1 = bulkm + 4. /3 * shearm;
    double a2 = bulkm - 2. /3 * shearm;

#ifdef THREED
    s[0] += a1 * de[0] + a2 * (de[1] + de[2]);
    s[1] += a1 * de[1] + a2 * (de[0] + de[2]);
    s[2] += a1 * de[2] + a2 * (de[0] + de[1]);
    s[3] += 2 * shearm * de[3];
    s[4] += 2 * shearm * de[4];
    s[5] += 2 * shearm * de[5];
#else
    s[0] += a1 * de[0] + a2 * de[1];
    s[1] += a1 * de[1] + a2 * de[0];
    s[2] += 2 * shearm * de[2];
#endif
}


static void maxwell(double bulkm, double shearm, double viscosity, double dt,
                    double dv, const double* de, double* s)
{
    // non-dimensional parameter: dt/ relaxation time
    double tmp = 0.5 * dt * shearm / viscosity;
    double f1 = 1 - tmp;
    double f2 = 1 / (1  + tmp);

    // deviatoric strain (diaganol component)
    double ded[NDIMS];
#ifdef THREED
    double dev = (de[0] + de[1] + de[2]) / 3;
    ded[0] = de[0] - dev;
    ded[1] = de[1] - dev;
    ded[2] = de[2] - dev;
#else
    double dev = (de[0] + de[1]) / 2;
    ded[0] = de[0] - dev;
    ded[1] = de[1] - dev;
#endif

    // deviatoric stress (diaganol component)
#ifdef THREED
    double s0 = (s[0] + s[1] + s[2]) / 3;
    s[0] -= s0;
    s[1] -= s0;
    s[2] -= s0;
#else
    double s0 = (s[0] + s[1]) / 2;
    s[0] -= s0;
    s[1] -= s0;
#endif

    // istropic stress is elastic
    s0 += bulkm * dv;

    // convert back to total stress
    for (int i=0; i<NDIMS; ++i)
        s[i] = (s[i] * f1 + 2 * shearm * ded[i]) * f2 + s0;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = (s[i] * f1 + 2 * shearm * de[i]) * f2;
}


void update_stress(const Variables& var, double2d& stress,
                   double2d& strain, double_vec& plstrain)
{
    const int rheol_type = var.mat->rheol_type;

    for (int e=0; e<var.nelem; ++e) {
        double bulkm, shearm;
        double viscosity;

        if (rheol_type & MatProps::rh_elastic) {
            bulkm = var.mat->bulkm(e);
            shearm = var.mat->shearm(e);
        }
        if (rheol_type & MatProps::rh_viscous) {
            viscosity = var.mat->visc(e);
            if (rheol_type == MatProps::rh_viscous)
                bulkm = var.mat->bulkm(e);
        }

        // TODO: strain correction

        // strain increment
        double de[NSTR];
        for (int i=0; i<NSTR; ++i) {
            de[i] = (*var.strain_rate)[e][i] * var.dt;
            strain[e][i] += de[i];
        }

        // stress of this element
        double* s = &stress[e][0];

        double dv;
        switch (rheol_type) {
        case MatProps::rh_elastic:
            elastic(bulkm, shearm, de, s);
            break;
        case MatProps::rh_viscous:
            std::cerr << "Error: pure viscous rheology not implemented.\n";
            std::exit(1);
            break;
        case MatProps::rh_maxwell:
            dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
            maxwell(bulkm, shearm, viscosity, var.dt, dv, de, s);
            break;
        default:
            std::cerr << "Error: unknown rheology type: " << rheol_type << "\n";
            std::exit(1);
            break;
        }
        std::cout << "stress " << e << ": ";
        print(std::cout, s, NSTR);
        std::cout << '\n';
    }
}
