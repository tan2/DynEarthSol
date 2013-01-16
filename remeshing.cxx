#include "iostream"

#include "constants.hpp"
#include "parameters.hpp"

#include "geometry.hpp"
#include "remeshing.hpp"


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
}
