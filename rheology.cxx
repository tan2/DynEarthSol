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

#ifdef THREED
    double dev = (de[0] + de[1] + de[2]) / 3;
    double s0 = (s[0] + s[1] + s[2]) / 3;
#else
    double dev = (de[0] + de[1]) / 2;
    double s0 = (s[0] + s[1]) / 2;
#endif

    // istropic stress is elastic
    s0 += bulkm * dv;

    // convert back to total stress
    for (int i=0; i<NDIMS; ++i)
        s[i] = ((s[i] - s0) * f1 + 2 * shearm * (de[i] - dev)) * f2 + s0;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = (s[i] * f1 + 2 * shearm * de[i]) * f2;
}


static void viscous(double bulkm, double viscosity, double total_dv,
                    const double* edot, double* s)
{
    /* Viscous Model + incompressibility enforced by bulk modulus */

#ifdef THREED
    double dev = (edot[0] + edot[1] + edot[2]) / 3;
#else
    double dev = (edot[0] + edot[1]) / 2;
#endif

    for (int i=0; i<NDIMS; ++i)
        s[i] = 2 * viscosity * (edot[i] - dev) + bulkm * total_dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = 2 * viscosity * edot[i];
}


static void elasto_plastic(double bulkm, double shearm,
                           double amc, double anphi, double anpsi,
                           double hardn, double ten_max,
                           const double* de, double& dpls, double* s)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion) */

    std::cerr << "Error: elastoc-plastic rheology not implemented.\n";
    std::exit(1);

}


void update_stress(const Variables& var, double2d& stress,
                   double2d& strain, double_vec& plstrain)
{
    const int rheol_type = var.mat->rheol_type;

    for (int e=0; e<var.nelem; ++e) {
        // stress, strain and strain_rate of this element
        double* s = &stress[e][0];
        double* es = &strain[e][0];
        const double *edot = &(*var.strain_rate)[e][0];

        // TODO: strain correction

        // strain increment
        double de[NSTR];
        for (int i=0; i<NSTR; ++i) {
            de[i] = es[i] * var.dt;
            es[i] += de[i];
        }

        switch (rheol_type) {
        case MatProps::rh_elastic:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                elastic(bulkm, shearm, de, s);
            }
            break;
        case MatProps::rh_viscous:
            {
                double bulkm = var.mat->bulkm(e);
                double viscosity = var.mat->visc(e);
#ifdef THREED
                double total_dv = es[0] + es[1] + es[2];
#else
                double total_dv = es[0] + es[1];
#endif
                viscous(bulkm, viscosity, total_dv, edot, s);
            }
            break;
        case MatProps::rh_maxwell:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double viscosity = var.mat->visc(e);
                double dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
                maxwell(bulkm, shearm, viscosity, var.dt, dv, de, s);
            }
            break;
        case MatProps::rh_ep:
            {
                double dpls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                               de, dpls, s);
                plstrain[e] += dpls;
            }
            break;
        case MatProps::rh_evp:
            {
                double dpls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double viscosity = var.mat->visc(e);
                double dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
                // stress due maxwell rheology
                double sv[NSTR];
                for (int i=0; i<NSTR; ++i) sv[i] = s[i];
                maxwell(bulkm, shearm, viscosity, var.dt, dv, de, sv);
                double svII = second_invariant2(sv);

                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                // stress due elasto-plastic rheology
                double sp[NSTR];
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                               de, dpls, sp);
                double spII = second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII)
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                else
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                plstrain[e] += dpls;
            }
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
