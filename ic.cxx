#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "ic.hpp"

namespace {

    class Zone
    {
    public:
        virtual ~Zone() {};
        virtual bool contains(const double x[NDIMS]) const = 0;
    };

    class Empty_zone : public Zone
    {
    public:
        bool contains(const double x[NDIMS]) const {return false;}
    };


    class Planar_zone : public Zone
    {
    private:
        const double az, incl;
        const double halfwidth; // in meter
#ifdef THREED
        const double ymin, ymax; // in meter
#endif
        const double zmin, zmax; // in meter
        const double *x0;

    public:
        Planar_zone(const double center[NDIMS], double azimuth, double inclination, double halfwidth_,
#ifdef THREED
                    double ymin_, double ymax_,
#endif
                    double zmin_, double zmax_) :
            az(std::tan(azimuth * DEG2RAD)), incl(1/std::tan(inclination * DEG2RAD)), halfwidth(halfwidth_),
#ifdef THREED
            ymin(ymin_), ymax(ymax_),
#endif
            zmin(zmin_), zmax(zmax_),
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {}

        bool contains(const double x[NDIMS]) const
        {
            // Is x within halfwidth distance to a plane cutting through x0?
            return (x[NDIMS-1] > zmin &&
                    x[NDIMS-1] < zmax &&
#ifdef THREED
                    x[1] > ymin &&
                    x[1] < ymax &&
#endif
                    std::fabs( (x[0] - x0[0])
#ifdef THREED
                               - az * (x[1] - x0[1])
#endif
                               + incl * (x[NDIMS-1] - x0[NDIMS-1]) ) < halfwidth );
        }
    };


    class Ellipsoidal_zone : public Zone
    {
    private:
        const double *x0;
        double semi_axis2[NDIMS];

    public:
        Ellipsoidal_zone(const double center[NDIMS], const double semi_axis[NDIMS]) :
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {
            for(int i=0; i<NDIMS; i++)
                semi_axis2[i] =  semi_axis[i] * semi_axis[i];
        }

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])/semi_axis2[0]
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])/semi_axis2[1]
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])/semi_axis2[NDIMS-1] < 1 );
        }
    };

} // anonymous namespace


double get_prem_pressure(double depth)
{
    // reference pressure profile from isotropic PREM model
    const int nlayers = 46;
    const double ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                                 60e3,   80e3,   115e3,  150e3,  185e3,
                                 220e3,  265e3,  310e3,  355e3,  400e3,
                                 450e3,  500e3,  550e3,  600e3,  635e3,
                                 670e3,  721e3,  771e3,  871e3,  971e3,
                                 1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                                 1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                                 2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                                 2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                                 2891e3 };

    // pressure in PREM table is given in kilobar, converted to 10^8 Pa
    const double ref_pressure[] = { 0e8,      0.3e8,    3.3e8,    6.0e8,    11.2e8,
                                    17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                                    71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                                    152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                                    238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                                    418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                                    655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                                    905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                                    1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                                    1357.5e8 };

    // PREM model doesn't work if depth is above sea level, always returns 0 pressure
    if (depth <= 0) return 0;

    int n;
    for (n=1; n<nlayers; n++) {
        if (depth <= ref_depth[n]) break;
    }

    // linear interpolation
    double pressure = ref_pressure[n-1] + (ref_pressure[n] - ref_pressure[n-1]) *
        (depth - ref_depth[n-1]) / (ref_depth[n] - ref_depth[n-1]);

    return pressure;
}


double get_prem_pressure_modified(double depth)
{
    // reference pressure profile from isotropic PREM model, modified for
    // average continental crust (density 2800 kg/m^3, thickness 24.4 km)
    const int nlayers = 46;
    const double ref_depth[] = { 0e3,    3e3,    15e3,   24.4e3, 40e3,
                                 60e3,   80e3,   115e3,  150e3,  185e3,
                                 220e3,  265e3,  310e3,  355e3,  400e3,
                                 450e3,  500e3,  550e3,  600e3,  635e3,
                                 670e3,  721e3,  771e3,  871e3,  971e3,
                                 1071e3, 1171e3, 1271e3, 1371e3, 1471e3,
                                 1571e3, 1671e3, 1771e3, 1871e3, 1971e3,
                                 2071e3, 2171e3, 2271e3, 2371e3, 2471e3,
                                 2571e3, 2671e3, 2741e3, 2771e3, 2871e3,
                                 2891e3 };

    // pressure in PREM table is given in kilobar, converted to 10^8 Pa
    const double ref_pressure[] = { 0e8,      0.82e8,    4.1e8,    6.7e8,    11.2e8,
                                    17.8e8,   24.5e8,   36.1e8,   47.8e8,   59.4e8,
                                    71.1e8,   86.4e8,   102.0e8,  117.7e8,  133.5e8,
                                    152.2e8,  171.3e8,  190.7e8,  210.4e8,  224.3e8,
                                    238.3e8,  260.7e8,  282.9e8,  327.6e8,  372.8e8,
                                    418.6e8,  464.8e8,  511.6e8,  558.9e8,  606.8e8,
                                    655.2e8,  704.1e8,  753.5e8,  803.6e8,  854.3e8,
                                    905.6e8,  957.6e8,  1010.3e8, 1063.8e8, 1118.2e8,
                                    1173.4e8, 1229.7e8, 1269.7e8, 1287.0e8, 1345.6e8,
                                    1357.5e8 };

    // PREM model doesn't work if depth is above sea level, always returns 0 pressure
    if (depth <= 0) return 0;

    int n;
    for (n=1; n<nlayers; n++) {
        if (depth <= ref_depth[n]) break;
    }

    // linear interpolation
    double pressure = ref_pressure[n-1] + (ref_pressure[n] - ref_pressure[n-1]) *
        (depth - ref_depth[n-1]) / (ref_depth[n] - ref_depth[n-1]);

    return pressure;
}


void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, tensor_t &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
            double p;
            if (param.control.ref_pressure_option == 1)
                p = get_prem_pressure(-zcenter);
            else if (param.control.ref_pressure_option == 2)
                p = get_prem_pressure_modified(-zcenter);
            for (int i=0; i<NDIMS; ++i) {
                stress[e][i] = -p;
                strain[e][i] = -p / ks / NDIMS;
            }
        }
        else if (param.control.ref_pressure_option == 0) {
            for (int i=0; i<NDIMS; ++i) {
                stress[e][i] = rho * param.control.gravity * zcenter;
                strain[e][i] = rho * param.control.gravity * zcenter / ks / NDIMS;
            }
        }
    }

    switch (param.control.ref_pressure_option) {
    case 0:
        compensation_pressure = rho * param.control.gravity * param.mesh.zlength;
        break;
    case 1:
        compensation_pressure = get_prem_pressure(param.mesh.zlength);
        break;
    case 2:
        compensation_pressure = get_prem_pressure_modified(param.mesh.zlength);
        break;
    default:
        std::cerr << "Error: unknown option for control.ref_pressure_option: " << param.control.ref_pressure_option << '\n';
        std::exit(1);
    }
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    Zone *weakzone;

    // TODO: adding different types of weak zone
    double plane_center[NDIMS]; // this variable must outlive weakzone
    switch (param.ic.weakzone_option) {
    case 0:
        weakzone = new Empty_zone();
        break;
    case 1:
        // a planar weak zone, cut through top center
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Planar_zone(plane_center,
                                   param.ic.weakzone_azimuth,
                                   param.ic.weakzone_inclination,
                                   param.ic.weakzone_halfwidth * param.mesh.resolution,
#ifdef THREED
                                   param.ic.weakzone_y_min * param.mesh.ylength,
                                   param.ic.weakzone_y_max * param.mesh.ylength,
#endif
                                   -param.ic.weakzone_depth_max * param.mesh.zlength,
                                   -param.ic.weakzone_depth_min * param.mesh.zlength);
        break;
    case 2:
        // a ellipsoidal weak zone
        double semi_axis[NDIMS];
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
        semi_axis[0] = param.ic.weakzone_xsemi_axis;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
        semi_axis[1] = param.ic.weakzone_ysemi_axis;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        semi_axis[NDIMS-1] = param.ic.weakzone_zsemi_axis;
        weakzone = new Ellipsoidal_zone(plane_center, semi_axis);
        break;
    default:
        std::cerr << "Error: unknown weakzone_option: " << param.ic.weakzone_option << '\n';
        std::exit(1);
    }

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[NDIMS] = {0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }

        if (weakzone->contains(center))
            plstrain[e] = param.ic.weakzone_plstrain;
    }

    delete weakzone;
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    const double age = param.ic.oceanic_plate_age_in_yr * YEAR2SEC;
    const MatProps &mat = *var.mat;
    const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element

    for (int i=0; i<var.nnode; ++i) {
        double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * age);
        temperature[i] = param.bc.surface_temperature +
            (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
    }
}


