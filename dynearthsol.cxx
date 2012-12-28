#include <iostream>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static void allocate_variables(Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    var.dvoldt = new double_vec(n);
    var.edvoldt = new double_vec(e);

    var.temperature = new double_vec(n);
    var.plstrain = new double_vec(e);
    var.tmp0 = new double_vec(std::max(n,e));

    var.vel = new arrayd2(n, 0);
    var.force = new arrayd2(n, 0);

    var.strain_rate = new tensord2(e, 0);
    var.strain = new tensord2(e, 0);
    var.stress = new tensord2(e, 0);

    var.shpdx = new shapefn(e);
    if (NDIMS == 3) var.shpdy = new shapefn(e);
    var.shpdz = new shapefn(e);
}


static void create_matprops(const Param &par, Variables &var)
{
    var.mat = new MatProps(par, var);
}


void compute_dvoldt(const Variables &var, double_vec &tmp,
                    double_vec &dvoldt, double_vec &edvoldt)
{
    /* dvoldt is the volumetric strain rate */
    /* edvoldt is the averaged dvoldt on the element */
    const double_vec& volume = *var.volume;
    const double_vec& volume_n = *var.volume_n;
    std::fill_n(tmp.begin(), var.nnode, 0);

    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                  \
            shared(egroup, var, tmp, volume)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = (*var.connectivity)[e];
        const double* strain_rate = (*var.strain_rate)[e];
        // TODO: try another definition:
        // dj = (volume[e] - volume_old[e]) / volume_old[e] / dt
        double dj = trace(strain_rate);
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            tmp[n] += dj * volume[e];
        }
    }
        } // end of ee
    }


    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, tmp, volume_n)
    for (int n=0; n<var.nnode; ++n)
         dvoldt[n] = tmp[n] / volume_n[n];

    #pragma omp parallel for default(none)      \
        shared(var, dvoldt, edvoldt)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double dj = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            int n = conn[i];
            dj += dvoldt[n];
        }
        edvoldt[e] = dj / NODES_PER_ELEM;
    }

    // std::cout << "dvoldt:\n";
    // print(std::cout, dvoldt);
    // std::cout << "\n";
    // std::cout << "edvoldt:\n";
    // print(std::cout, edvoldt);
    // std::cout << "\n";
}


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

    // pressure in PREM table is given in kilobar, convert to 10^8 Pa
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

    int n;
    for (n=1; n<nlayers; n++) {
        if (depth <= ref_depth[n]) break;
    }

    // linear interpolation
    double pressure = ref_pressure[n-1] + depth *
        (ref_pressure[n] - ref_pressure[n-1]) / (ref_depth[n] - ref_depth[n-1]);

    return std::max(pressure, 0.0);
}


void initial_stress_state(const Param &param, const Variables &var,
                          tensord2 &stress, tensord2 &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    int earthlike_reference_pressure = 0;
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        zcenter /= NODES_PER_ELEM;

        if (param.control.ref_pressure_option == 1) {
            ks = var.mat->bulkm(e);
            double p = get_prem_pressure(-zcenter);
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
    default:
        std::cerr << "Error: unknown option for control.ref_pressure_option: " << param.control.ref_pressure_option << '\n';
        std::exit(1);
    }
}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[3] = {0,0,0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }
        // TODO: adding different types of weak zone
        const double d = param.mesh.resolution;
        if (std::fabs(center[0] - param.mesh.xlength * 0.5) < 2*d &&
            std::fabs(center[NDIMS-1] + param.mesh.zlength) < 1.5*d)
            plstrain[e] = 0.1;
    }
}


void initial_temperature(const Param &param, const Variables &var, double_vec &temperature)
{
    const double oceanic_plate_age = 60e6 * YEAR2SEC;
    const double diffusivity = 1e-6;

    for (int i=0; i<var.nnode; ++i) {
        double w = -(*var.coord)[i][NDIMS-1] / std::sqrt(4 * diffusivity * oceanic_plate_age);
        temperature[i] = param.bc.surface_temperature +
            (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
    }
}


static bool is_on_boundary(const Variables &var, int node)
{
    int flag = (*var.bcflag)[node];
    return flag & (BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | BOUNDZ0 | BOUNDZ1);
}


static double find_max_vbc(const BC &bc)
{
    double max_vbc_val = 0;
    if (bc.vbc_x0 == 1 || bc.vbc_x0 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x0));
    if (bc.vbc_x1 == 1 || bc.vbc_x1 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_x1));
    if (bc.vbc_y0 == 1 || bc.vbc_y0 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y0));
    if (bc.vbc_y1 == 1 || bc.vbc_y1 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_y1));
    if (bc.vbc_z0 == 1 || bc.vbc_z0 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z0));
    if (bc.vbc_z1 == 1 || bc.vbc_z1 == 3)
        max_vbc_val = std::max(max_vbc_val, std::fabs(bc.vbc_val_z1));

    return max_vbc_val;
}


void apply_vbcs(const Param &param, const Variables &var, arrayd2 &vel)
{
    // meaning of vbc flags (odd: free; even: fixed) --
    // 0: all components free
    // 1: normal component fixed, shear components free
    // 2: normal component free, shear components fixed at 0
    // 3: normal component fixed, shear components fixed at 0
    // 4: normal component free, shear component (not z) fixed, only in 3D
    // 5: normal component fixed at 0, shear component (not z) fixed, only in 3D

    const BC &bc = param.bc;

    // diverging x-boundary
    #pragma omp parallel for default(none) \
        shared(bc, var, vel)
    for (int i=0; i<var.nnode; ++i) {

        // fast path: skip nodes not on boundary
        if (! is_on_boundary(var, i)) continue;

        int flag = (*var.bcflag)[i];
        double *v = vel[i];

        // X
        if (flag & BOUNDX0) {
            switch (bc.vbc_x0) {
            case 0:
                break;
            case 1:
                v[0] = bc.vbc_val_x0;
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
                v[0] = bc.vbc_val_x0;
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc.vbc_val_x0;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc.vbc_val_x0;
                v[2] = 0;
                break;
#endif
            }
        }
        else if (flag & BOUNDX1) {
            switch (bc.vbc_x1) {
            case 0:
                break;
            case 1:
                v[0] = bc.vbc_val_x1;
                break;
            case 2:
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
            case 3:
                v[0] = bc.vbc_val_x1;
                v[1] = 0;
#ifdef THREED
                v[2] = 0;
#endif
                break;
#ifdef THREED
            case 4:
                v[1] = bc.vbc_val_x1;
                v[2] = 0;
                break;
            case 5:
                v[0] = 0;
                v[1] = bc.vbc_val_x1;
                v[2] = 0;
                break;
#endif
            }
        }
#ifdef THREED
        // Y
        if (flag & BOUNDY0) {
            switch (bc.vbc_y0) {
            case 0:
                break;
            case 1:
                v[1] = bc.vbc_val_y0;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[1] = bc.vbc_val_y0;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            case 5:
                v[1] = 0;
                v[0] = bc.vbc_val_y0;
                v[2] = 0;
                break;
            }
        }
        else if (flag & BOUNDY1) {
            switch (bc.vbc_y1) {
            case 0:
                break;
            case 1:
                v[1] = bc.vbc_val_y1;
                break;
            case 2:
                v[0] = 0;
                v[2] = 0;
                break;
            case 3:
                v[1] = bc.vbc_val_y1;
                v[1] = 0;
                v[2] = 0;
                break;
            case 4:
                v[0] = bc.vbc_val_y1;
                v[2] = 0;
                break;
            case 5:
                v[1] = 0;
                v[0] = bc.vbc_val_y1;
                v[2] = 0;
                break;
            }
        }
#endif
        // Z

        // fast path: vz is usually free in the models
        if (bc.vbc_z0==0 && bc.vbc_z1==0) continue;

        if (flag & BOUNDZ0) {
            switch (bc.vbc_z0) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc.vbc_val_z0;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc.vbc_val_z0;
                break;
            case 4:
                v[0] = bc.vbc_val_z0;
#ifdef THREED
                v[1] = bc.vbc_val_z0;
#endif
                break;
            case 5:
                v[0] = bc.vbc_val_z0;
#ifdef THREED
                v[1] = bc.vbc_val_z0;
#endif
                v[NDIMS-1] = 0;
                break;
            }
        }
        else if (flag & BOUNDZ1) {
            switch (bc.vbc_z1) {
            case 0:
                break;
            case 1:
                v[NDIMS-1] = bc.vbc_val_z1;
                break;
            case 2:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                break;
            case 3:
                v[0] = 0;
#ifdef THREED
                v[1] = 0;
#endif
                v[NDIMS-1] = bc.vbc_val_z1;
                break;
            case 4:
                v[0] = bc.vbc_val_z1;
#ifdef THREED
                v[1] = bc.vbc_val_z1;
#endif
                break;
            case 5:
                v[0] = bc.vbc_val_z1;
#ifdef THREED
                v[1] = bc.vbc_val_z1;
#endif
                v[NDIMS-1] = 0;
                break;
            }
        }
    }
}


void init(const Param& param, Variables& var)
{
    create_new_mesh(param, var);
    allocate_variables(var);
    create_matprops(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    *var.volume_old = *var.volume;
    compute_dvoldt(var, *var.tmp0, *var.dvoldt, *var.edvoldt);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    initial_stress_state(param, var, *var.stress, *var.strain, var.compensation_pressure);
    initial_weak_zone(param, var, *var.plstrain);
    initial_temperature(param, var, *var.temperature);
    apply_vbcs(param, var, *var.vel);
}


void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{
    // diffusion matrix
    double D[NODES_PER_ELEM][NODES_PER_ELEM];

    tdot.assign(var.nnode, 0);
    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                                  \
            shared(egroup, var, param, temperature, tdot) private(D)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = (*var.connectivity)[e];
        double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volumn
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int j=0; j<NODES_PER_ELEM; ++j) {
#ifdef THREED
                D[i][j] = (shpdx[i] * shpdx[j] +
                           shpdy[i] * shpdy[j] +
                           shpdz[i] * shpdz[j]);
#else
                D[i][j] = (shpdx[i] * shpdx[j] +
                           shpdz[i] * shpdz[j]);
#endif
            }
        }
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double diffusion = 0;
            for (int j=0; j<NODES_PER_ELEM; ++j)
                diffusion += D[i][j] * temperature[conn[j]];

            tdot[conn[i]] += diffusion * kv;
        }
    }
        }
    }

     #pragma omp parallel for default(none)      \
         shared(var, param, tdot, temperature)
     for (int n=0; n<var.nnode; ++n) {
        if ((*var.bcflag)[n] & BOUNDZ1)
            temperature[n] = param.bc.surface_temperature;
        else
            temperature[n] -= tdot[n] * var.dt / (*var.tmass)[n];
    }
}


void update_strain_rate(const Variables& var, tensord2& strain_rate)
{
    double *v[NODES_PER_ELEM];

    #pragma omp parallel for default(none) \
        shared(var, strain_rate) private(v)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.strain_rate)[e];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = (*var.vel)[conn[i]];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
        // YY component
        n = 1;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][1] * shpdy[i];
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][NDIMS-1] * shpdz[i];

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] + v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] + v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] + v[i][2] * shpdy[i]);
#endif
    }
}


static void apply_stress_bcs(const Param& param, const Variables& var, arrayd2& force)
{
    // TODO: add general stress (Neumann) bcs

    // TODO: add water loading from the surface boundary

    // Wrinkler foundation for the bottom boundary
    if (param.bc.wrinkler_foundation && param.control.gravity != 0) {
        const int bottom_bdry = bdry_order.find(BOUNDZ0)->second;
        const auto& bottom = var.bfacets[bottom_bdry];
        const auto& coord = *var.coord;
        // loops over all bottom facets
        for (std::size_t i=0; i<bottom.size(); ++i) {
            // this facet belongs to element e
            int e = bottom[i].first;
            // this facet is the f-th facet of e
            int f = bottom[i].second;
            const int *conn = (*var.connectivity)[e];

            int n0 = conn[NODE_OF_FACET[f][0]];
            int n1 = conn[NODE_OF_FACET[f][1]];
#ifdef THREED
            int n2 = conn[NODE_OF_FACET[f][2]];

            // vectors: n0-n1 and n0-n2
            double v01[NDIMS], v02[NDIMS];
            for (int i=0; i<NDIMS; ++i) {
                v01[i] = coord[n1][i] - coord[n0][i];
                v02[i] = coord[n2][i] - coord[n0][i];
            }

            // the outward normal vector is parallel to the cross product
            // of v01 and v02, with magnitude = area of the triangle
            double normal[NDIMS];
            normal[0] = (v01[1] * v02[2] - v01[2] * v02[1]) / 2;
            normal[1] = (v01[2] * v02[0] - v01[0] * v02[2]) / 2;
            normal[2] = (v01[0] * v02[1] - v01[1] * v02[0]) / 2;

            double zcenter = (coord[n0][2] + coord[n1][2] + coord[n2][2]) / NODES_PER_FACET;
#else
            // the edge vector
            double v01[NDIMS];
            for (int i=0; i<NDIMS; ++i) {
                v01[i] = coord[n1][i] - coord[n0][i];
            }

            // the normal vector to the edge, pointing outward
            double normal[NDIMS];
            normal[0] = v01[1];
            normal[1] = -v01[0];

            double zcenter = (coord[n0][1] + coord[n1][1]) / NODES_PER_FACET;
#endif
            double dz = zcenter - (-param.mesh.zlength);
            double p = var.compensation_pressure -
                (var.mat->rho(e) + param.bc.wrinkler_delta_rho) * param.control.gravity * dz;

            // bottom support - Archimed force (normal to the surface)
            for (int i=0; i<NDIMS; ++i) {
                force[n0][i] -= p * normal[i] / NODES_PER_FACET;
                force[n1][i] -= p * normal[i] / NODES_PER_FACET;
#ifdef THREED
                force[n2][i] -= p * normal[i] / NODES_PER_FACET;
#endif
            }
        }
    }
}


static void apply_damping(const Param& param, const Variables& var, arrayd2& force)
{
    // flatten 2d arrays to simplify indexing
    double* ff = force.data();
    const double* v = var.vel->data();
    const double small_vel = 1e-13;
    #pragma omp parallel for default(none)          \
        shared(var, param, ff, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        if (std::fabs(v[i]) > small_vel) {
            ff[i] -= param.control.damping_factor * std::copysign(ff[i], v[i]);
        }
    }
}


void update_force(const Param& param, const Variables& var, arrayd2& force)
{
    std::fill_n(force.data(), var.nnode*NDIMS, 0);

    for (auto egroup=var.egroups->begin(); egroup!=var.egroups->end(); egroup++) {
        #pragma omp parallel for default(none)                  \
            shared(egroup, param, var, force)
        for (std::size_t ee=0; ee<egroup->size(); ++ee) {
            int e = (*egroup)[ee];
    {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
#endif
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.stress)[e];
        double vol = (*var.volume)[e];

        double buoy = 0;
        if (param.control.gravity != 0)
            buoy = var.mat->rho(e) * param.control.gravity / NODES_PER_ELEM;

        for (int i=0; i<NODES_PER_ELEM; ++i) {
            double *f = force[conn[i]];
#ifdef THREED
            f[0] -= (s[0]*shpdx[i] + s[3]*shpdy[i] + s[4]*shpdz[i]) * vol;
            f[1] -= (s[3]*shpdx[i] + s[1]*shpdy[i] + s[5]*shpdz[i]) * vol;
            f[2] -= (s[4]*shpdx[i] + s[5]*shpdy[i] + s[2]*shpdz[i] + buoy) * vol;
#else
            f[0] -= (s[0]*shpdx[i] + s[2]*shpdz[i]) * vol;
            f[1] -= (s[2]*shpdx[i] + s[1]*shpdz[i] + buoy) * vol;
#endif
        }
    }
        } // end of ee
    }

    apply_stress_bcs(param, var, force);

    if (param.control.damping_factor != 0) {
        apply_damping(param, var, force);
    }
}


void rotate_stress() {};


void update_velocity(const Variables& var, arrayd2& vel)
{
    const double* m = &(*var.mass)[0];
    // flatten 2d arrays to simplify indexing
    const double* f = var.force->data();
    double* v = vel.data();
    #pragma omp parallel for default(none) \
        shared(var, m, f, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        int n = i / NDIMS;
        v[i] += var.dt * f[i] / m[n];
    }
}


static void update_coordinate(const Variables& var, arrayd2& coord)
{
    double* x = var.coord->data();
    const double* v = var.vel->data();

    #pragma omp parallel for default(none) \
        shared(var, x, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        x[i] += v[i] * var.dt;
    }

    // surface_processes()
}


void update_mesh(const Param& param, Variables& var)
{
    update_coordinate(var, *var.coord);

    var.volume->swap(*var.volume_old);
    compute_volume(*var.coord, *var.connectivity, *var.egroups, *var.volume, *var.volume_n);
    compute_dvoldt(var, *var.tmp0, *var.dvoldt, *var.edvoldt);
    compute_mass(param, *var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.mass, *var.tmass);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);
}


int main(int argc, const char* argv[])
{
    double start_time = 0;
#ifdef USE_OMP
    start_time = omp_get_wtime();
#endif

    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        std::cout << "       " << argv[0] << " -h or --help\n";
        return -1;
    }

    Param param;
    void get_input_parameters(const char*, Param&);
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    var.time = 0;
    var.steps = 0;
    var.frame = 0;

    var.max_vbc_val = find_max_vbc(param.bc);

    if (! param.sim.is_restarting) {
        init(param, var);
        output(param, var, start_time);
        var.frame ++;
    }
    else {
        restart();
        var.frame ++;
    }

    var.dt = compute_dt(param, var);

    do {
        var.steps ++;
        var.time += var.dt;

        update_temperature(param, var, *var.temperature, *var.tmp0);
        update_strain_rate(var, *var.strain_rate);
        update_stress(var, *var.stress, *var.strain, *var.plstrain, *var.strain_rate);
        update_force(param, var, *var.force);
        update_velocity(var, *var.vel);
        apply_vbcs(param, var, *var.vel);
        update_mesh(param, var);
        // dt computation is expensive, and dt only changes slowly
        // don't have to do it every time step
        if (var.steps % 10 == 0) var.dt = compute_dt(param, var);
        rotate_stress();

        if ( (var.steps == var.frame * param.sim.output_step_interval) ||
             (var.time > var.frame * param.sim.output_time_interval_in_yr * YEAR2SEC) ) {
            output(param, var, start_time);
            var.frame ++;
        }

    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    return 0;
}
