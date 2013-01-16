#include "iostream"

#include "constants.hpp"
#include "parameters.hpp"

#include "geometry.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "remeshing.hpp"


static void reallocate_variables(Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    delete var.volume;
    delete var.volume_old;
    delete var.volume_n;
    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    delete var.mass;
    delete var.tmass;
    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    delete var.dvoldt;
    delete var.edvoldt;
    var.dvoldt = new double_vec(n);
    var.edvoldt = new double_vec(e);

    delete var.temperature;
    delete var.plstrain;
    delete var.tmp0;
    var.temperature = new double_vec(n);
    var.plstrain = new double_vec(e);
    var.tmp0 = new double_vec(std::max(n,e));

    delete var.vel;
    delete var.force;
    var.vel = new array_t(n, 0);
    var.force = new array_t(n, 0);

    delete var.strain_rate;
    delete var.strain;
    delete var.stress;
    var.strain_rate = new tensor_t(e, 0);
    var.strain = new tensor_t(e, 0);
    var.stress = new tensor_t(e, 0);

    delete var.shpdx;
    delete var.shpdz;
    var.shpdx = new shapefn(e);
    if (NDIMS == 3) {
        delete var.shpdy;
        var.shpdy = new shapefn(e);
    }
    var.shpdz = new shapefn(e);
}


bool bad_mesh_quality(const Param &param, const Variables &var)
{
    int worst_elem;
    double q = worst_elem_quality(*var.coord, *var.connectivity,
                                  *var.volume, worst_elem);
#ifdef THREED
    // normalizing q so that its magnitude is about the same in 2D and 3D
    q = std::pow(q, 1.0/3);
#endif
    std::cout << "Worst mesh quality = " << q << " at element #" << worst_elem << ".\n";
    if (q < param.mesh.min_quality) {
        return 1;
    }
    return 0;
}


void remesh(const Param &param, Variables &var)
{
    std::cout << "  Remeshing starts...\n";

    // setting up barycentric transformation

    // saving the mesh to .poly file

    // creating a "copy" of mesh pointer so that they are not deleted
    array_t old_coord;
    conn_t old_connectivity;
    segment_t old_segment;
    segflag_t old_segflag;
    old_coord.steal_ref(*var.coord);
    old_connectivity.steal_ref(*var.connectivity);
    old_segment.steal_ref(*var.segment);
    old_segflag.steal_ref(*var.segflag);

    delete var.coord;
    delete var.connectivity;
    delete var.segment;
    delete var.segflag;

    // XXX: modifying boundary if necessary

    // deleting (non-boundary) nodes to avoid small elements

    // saving the modified mesh to .poly file

    // new mesh
    int npoints = var.nnode;
    double *points = old_coord.data();
    int n_init_segments = var.nseg;
    int *init_segments = old_segment.data();
    int *init_segflags = old_segflag.data();

    // We don't want to refine large elements during remeshing,
    // so using the domain size as the max area
    double max_elem_size;
#ifdef THREED
    max_elem_size = param.mesh.xlength * param.mesh.ylength * param.mesh.zlength;
#else
    max_elem_size = param.mesh.xlength * param.mesh.zlength;
#endif
    double vertex_per_polygon = 3;
    points_to_mesh(param, var, npoints, points,
                   n_init_segments, init_segments, init_segflags,
                   max_elem_size, vertex_per_polygon);

    // interpolating fields

    // arrays of new mesh
    reallocate_variables(var);
    delete var.mat;
    var.mat = new MatProps(param, var);

    // updating other arrays

    free(old_coord.data());
    free(old_connectivity.data());
    free(old_segment.data());
    free(old_segflag.data());

    std::cout << "  Remeshing finished.\n";
}


