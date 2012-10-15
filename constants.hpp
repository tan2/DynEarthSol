#ifndef DYNEARTHSOL3D_CONSTANTS_HPP
#define DYNEARTHSOL3D_CONSTANTS_HPP

# include <cmath>

// # of spatial dimension, either 2 or 3
const int NDIMS = 2;

// triangles (3 node) in 2D and tetrahedra (4 nodes) in 3D
const int NODES_PER_ELEM = NDIMS + 1;

// Flags for boundary
const int BOUNDX0 = 1;
const int BOUNDX1 = 2;
const int BOUNDY0 = 3;
const int BOUNDY1 = 4;
const int BOUNDZ0 = 5;
const int BOUNDZ1 = 6;


const double YEAR2SEC = 365.2422 * 86400;
const double DEG2RAD = M_PI / 180;

#endif
