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
    std::cout << "Worst mesh quality = " << q << " at element #" << worst_elem << ".\n";
    if (q < 0.4) {
        return 1;
    }
    return 0;
}


void remesh(const Param &param, Variables &var)
{
}
