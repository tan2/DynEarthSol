#include <cmath>
#include <iostream>

#include "3x3-C/dsyevh3.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static inline
double trace(const double* s)
{
#ifdef THREED
    return s[0] + s[1] + s[2];
#else
    return s[0] + s[1];
#endif
}


static void principal_stresses3(const double* s, double p[3], double v[3][3])
{
    /* s is a flattened stress vector, with the components {XX, YY, ZZ, XY, XZ, YZ}.
     * Returns the eigenvalues p and eignvectors v.
     * The eigenvalues are ordered such that p[0] <= p[1] <= p[2].
     */

    // unflatten s to a 3x3 tensor, only the upper part is needed.
    double a[3][3];
    a[0][0] = s[0];
    a[1][1] = s[1];
    a[2][2] = s[2];
    a[0][1] = s[3];
    a[0][2] = s[4];
    a[1][2] = s[5];

    dsyevh3(a, v, p);
}


static void principal_stresses2(const double* s, double syy, double p[3],
                                double& cos2t, double& sin2t, int& icase)
{
    /* 's' is a flattened stress vector, with the components {XX, ZZ, XZ}.
     * The YY component is passed in in 'syy', assumed to be one of eigenvalue.
     * Returns the eigenvalues 'p', and the direction cosine of the
     * eigenvectors in the X-Z plane.
     * The eigenvalues are ordered such that p[0] <= p[1] <= p[2].
     */

    // center and radius of Mohr circle
    double s0 = 0.5 * (s[0] + s[1]);
    double rad = second_invariant(s);

    // principal stresses in the X-Z plane
    double si = s0 - rad;
    double sii = s0 + rad;

    {
        // direction cosine and sine of 2*theta
        const double eps = 1e-15;
        double a = s[0] - s[1];
        double b = - 2 * rad; // always negative
        if (b < -eps) {
            cos2t = a / b;
            sin2t = 2 * s[2] / b;
        }
        else {
            cos2t = 1;
            sin2t = 0;
        }
    }

    // sort p.s.
    if (syy > sii) {
        // syy is minor p.s.
        icase = 3;
        p[0] = si;
        p[1] = sii;
        p[2] = syy;
    }
    else if (syy < si) {
        // syy is major p.s.
        icase = 2;
        p[0] = syy;
        p[1] = si;
        p[2] = sii;
    }
    else {
        // syy is intermediate
        icase = 1;
        p[0] = si;
        p[1] = syy;
        p[2] = sii;
    }
}


static void elastic(double bulkm, double shearm, const double* de, double* s)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + lambda * dev;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];
}


static void maxwell(double bulkm, double shearm, double viscosity, double dt,
                    double dv, const double* de, double* s)
{
    // non-dimensional parameter: dt/ relaxation time
    double tmp = 0.5 * dt * shearm / viscosity;
    double f1 = 1 - tmp;
    double f2 = 1 / (1  + tmp);

    double dev = trace(de) / NDIMS;
    double s0 = trace(s) / NDIMS;

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

    double dev = trace(edot) / NDIMS;

    for (int i=0; i<NDIMS; ++i)
        s[i] = 2 * viscosity * (edot[i] - dev) + bulkm * total_dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = 2 * viscosity * edot[i];
}


static void elasto_plastic(double bulkm, double shearm,
                           double amc, double anphi, double anpsi,
                           double hardn, double ten_max,
                           const double* de, double& depls, double* s)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion) */

    // elastic trial stress
    elastic(bulkm, shearm, de, s);
    depls = 0;

    //
    // transform to principal stress coordinate system
    //
    // eigenvalues (principal stresses)
    double p[3];
#ifndef THREED
    // Stress YY component, plain strain
    double syy = (bulkm - 2./3 * shearm) * (de[1] + de[2]);
    // In 2D, we only construct the eigenvectors from
    // cos(2*theta) and sin(2*theta) of Mohr circle
    double cos2t, sin2t;
    int icase;
    principal_stresses2(s, syy, p, cos2t, sin2t, icase);
#else
    // eigenvectors
    double v[3][3];
    principal_stresses3(s, p, v);
#endif

    // composite (shear and tensile) yield criterion
    double fs = p[0] - p[NDIMS-1] * anphi + amc;
    double ft = p[NDIMS-1] - ten_max;

    if (fs > 0 && ft < 0) {
        // no failure
        return;
    }

    // yield, shear or tensile?
    double pa = std::sqrt(1 + anphi*anphi) + anphi;
    double ps = ten_max * anphi - amc;
    double h = p[NDIMS-1] - ten_max + pa * (p[0] - ps);
    double a1 = bulkm + 4. / 3 * shearm;
    double a2 = bulkm - 2. / 3 * shearm;

    double alam;
    if (h < 0) {
        // shear failure
        alam = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + hardn);
        p[0] -= alam * (a1 - a2 * anpsi);
        p[1] -= alam * (a2 - a2 * anpsi);
        p[2] -= alam * (a2 - a1 * anpsi);

        // 2nd invariant of plastic strain
#ifdef THREED
        /* // plastic strain in the principle directions, depls2 is always 0
         * double depls1 = alam;
         * double depls3 = -alam * anpsi;
         * double deplsm = (depls1 + depls3) / 3;
         * depls = std::sqrt( (depls1-deplsm)*(depls1-deplsm) +
         *                   (-deplsm)*(-deplsm) +
         *                   (depls3-deplsm)*(depls3-deplsm) );
         */
        // the equations above can be reduce to:
        depls = std::fabs(alam) * std::sqrt((1 + anpsi + anpsi*anpsi) / 3);
#else
        // TODO: copied from flac, but it is not consistent with the 3D formula above.
        // check with Luc and Eunseo ...
        depls = 0.5 * std::fabs(alam + alam * anpsi);
#endif
    }
    else {
        // tensile failure
        // TODO: the final stress in the tensile failure case, copied from Snac,
        // but the formula is different from the formula in flac.
        // check with Luc and Eunseo ...
        alam = ft / a1;
        p[0] -= alam * a2;
        p[1] -= alam * a2;
        p[2] -= alam * a1;

        // 2nd invariant of plastic strain
#ifdef THREED
        depls = std::fabs(alam) / std::sqrt(3);
#else
        // TODO: flac doesn't define depls in the tensile failure case.
        // check with Luc and Eunseo ...
        depls = alam;
#endif
    }

    // rotate the principal stresses back to global axes
    {
#ifdef THREED
        double ss[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for(int m=0; m<3; m++) {
            for(int n=m; n<3; n++) {
                for(int k=0; k<3; k++) {
                    ss[m][n] += v[k][m] * v[k][n] * p[k];
                }
            }
        }
        s[0] = ss[0][0];
        s[1] = ss[1][1];
        s[2] = ss[2][2];
        s[3] = ss[0][1];
        s[4] = ss[0][2];
        s[5] = ss[1][2];
#else
        int n1, n2;
        switch (icase) {
        case 1:
            n1 = 0;
            n2 = 2;
            break;
        case 2:
            n1 = 1;
            n2 = 2;
            break;
        case 3:
            n1 = 0;
            n2 = 1;
            break;
        }
        double dc2 = (p[n1] - p[n2]) * cos2t;
        double dss = p[n1] + p[n2];
        s[0] = 0.5 * (dss + dc2);
        s[1] = 0.5 * (dss - dc2);
        s[2] = 0.5 * (p[n1] - p[n2]) * sin2t;
#endif
    }
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
            de[i] = edot[i] * var.dt;
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
                double total_dv = trace(es);
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
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                               de, depls, s);
                plstrain[e] += depls;
            }
            break;
        case MatProps::rh_evp:
            {
                double depls = 0;
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
                               de, depls, sp);
                double spII = second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII)
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                else
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                plstrain[e] += depls;
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
