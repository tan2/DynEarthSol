#ifndef DYNEARTHSOL3D_CONSTANTS_HPP
#define DYNEARTHSOL3D_CONSTANTS_HPP

# include <cmath>
# include <map>

// # of spatial dimension, either 2 or 3
#ifdef THREED
const int NDIMS = 3;
#else
const int NDIMS = 2;
#endif

// triangles (3 node) in 2D and tetrahedra (4 nodes) in 3D
const int NODES_PER_ELEM = NDIMS + 1;

// # of indep. components of a symmetric tensor
// 2D -> 3;  3D -> 6
const int NSTR = NDIMS * (NDIMS + 1) / 2;

// Flags for boundary
const int iboundx0 = 0;
const int iboundx1 = 1;
const int iboundy0 = 2;
const int iboundy1 = 3;
const int iboundz0 = 4;
const int iboundz1 = 5;

typedef unsigned int uint;
const uint BOUNDX0 = 1 << iboundx0;  //  1, left
const uint BOUNDX1 = 1 << iboundx1;  //  2, right
const uint BOUNDY0 = 1 << iboundy0;  //  4, front
const uint BOUNDY1 = 1 << iboundy1;  //  8, back
const uint BOUNDZ0 = 1 << iboundz0;  // 16, bottom
const uint BOUNDZ1 = 1 << iboundz1;  // 32, top

// # of facets (edges) per element, 3 for 2D, 4 for 3D
const int FACETS_PER_ELEM = NDIMS + 1;
// # of nodes per facets (edges), 2 for 2D, 3 for 3D
const int NODES_PER_FACET = NDIMS;

// local node # of a facet
#ifdef THREED
// the nodes are ordered counter-clockwise, if viewed from outside
const int NODE_OF_FACET[FACETS_PER_ELEM][NODES_PER_FACET] =
    {{1,2,3},
     {0,3,2},
     {0,1,3},
     {0,2,1}};
#else
// the nodes are ordered counter-clockwise, if viewed from above
const int NODE_OF_FACET[FACETS_PER_ELEM][NODES_PER_FACET] =
    {{1,2},
     {2,0},
     {0,1}};
#endif

const double YEAR2SEC = 365.2422 * 86400;
const double DEG2RAD = M_PI / 180;

#endif
