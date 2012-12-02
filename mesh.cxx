
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef THREED

#define TETLIBRARY
#include "tetgen/tetgen.h"
#undef TETLIBRARY

#else

#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS

#endif // THREED

#include "constants.hpp"
#include "parameters.hpp"
#include "sortindex.hpp"
#include "utils.hpp"
#include "mesh.hpp"


void triangulate_polygon
(double min_angle, double max_area,
 int npoints, int nsegments,
 double *points, int *segments, int *segflags,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags)
{
#ifndef THREED
    char options[255];
    triangulateio in, out;

    // add 'Q' for no output; add multiple 'V's for verbose output
    char verbosity[] = "V";
    std::sprintf(options, "%spq%fjza%f", verbosity, min_angle, max_area);
    //std::puts(options);

    in.pointlist = points;
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.numberofpoints = npoints;
    in.numberofpointattributes = 0;

    in.trianglelist = NULL;
    in.triangleattributelist = NULL;
    in.trianglearealist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 3;
    in.numberoftriangleattributes = 0;

    in.segmentlist = segments;
    in.segmentmarkerlist = segflags;
    in.numberofsegments = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.regionlist = NULL;
    in.numberofregions = 0;

    out.pointlist = NULL;
    out.pointattributelist = NULL;
    out.pointmarkerlist = NULL;
    out.trianglelist = NULL;
    out.triangleattributelist = NULL;
    out.neighborlist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.edgelist = NULL;
    out.edgemarkerlist = NULL;

    /*******************************/
    triangulate(options, &in, &out, NULL);
    /*******************************/

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;

    *ntriangles = out.numberoftriangles;
    *triangles = out.trianglelist;

    *noutsegments = out.numberofsegments;
    *outsegments = out.segmentlist;
    *outsegflags = out.segmentmarkerlist;

    trifree(out.pointmarkerlist);
#endif
}


void tetrahedralize_polyhedron
(double max_ratio, double min_dihedral_angle, double max_volume,
 int npoints, int nsegments,
 double *points, int *segments, int *segflags,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags)
{
#ifdef THREED
    //
    // Setting Tetgen options.
    //
    char options[255];
    double max_dihedral_angle = 180 - 3 * min_dihedral_angle;

    // add 'Q' for no output; add multiple 'V's for verbose output
    char verbosity[] = "V";
    std::sprintf(options, "%spzq%fqq%fqqq%fa%f", verbosity, max_ratio,
                 min_dihedral_angle, max_dihedral_angle, max_volume);
    //std::puts(options);

    //
    // Setting input arrays to tetgen
    //
    tetgenio in;
    in.pointlist = points;
    in.numberofpoints = npoints;

    tetgenio::polygon *polys = new tetgenio::polygon[nsegments];
    for (int i=0; i<nsegments; ++i) {
        polys[i].vertexlist = &segments[i*4];
        polys[i].numberofvertices = 4;
    }

    tetgenio::facet *fl = new tetgenio::facet[nsegments];
    for (int i=0; i<nsegments; ++i) {
        fl[i].polygonlist = &polys[i];
        fl[i].numberofpolygons = 1;
        fl[i].holelist = NULL;
        fl[i].numberofholes = 0;
    }

    in.facetlist = fl;
    in.facetmarkerlist = segflags;
    in.numberoffacets = nsegments;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.regionlist = NULL;
    in.numberofregions = 0;

    tetgenio out;
    /*******************************/
    tetrahedralize(options, &in, &out, NULL, NULL);
    /*******************************/

    // the destructor of tetgenio will free any non-NULL pointer
    // set in.pointers to NULL to prevent double-free
    in.pointlist = NULL;
    in.facetmarkerlist = NULL;
    in.facetlist = NULL;
    delete [] polys;
    delete [] fl;

    *noutpoints = out.numberofpoints;
    *outpoints = out.pointlist;
    out.pointlist = NULL;

    *ntriangles = out.numberoftetrahedra;
    *triangles = out.tetrahedronlist;
    out.tetrahedronlist = NULL;

    *noutsegments = out.numberoftrifaces;
    *outsegments = out.trifacelist;
    *outsegflags = out.trifacemarkerlist;
    out.trifacelist = NULL;
    out.trifacemarkerlist = NULL;

#endif
}


static void new_mesh_uniform_resolution(const Param& param, Variables& var)
{
    int npoints = 4 * (NDIMS - 1); // 2D:4;  3D:8
    double *points = new double[npoints*NDIMS];

    int n_init_segments = 2 * NDIMS; // 2D:4;  3D:6
    int n_segment_nodes = 2 * (NDIMS - 1); // 2D:2; 3D:4
    int *init_segments = new int[n_init_segments*n_segment_nodes];
    int *init_segflags = new int[n_init_segments];

    int nnode, nelem, nseg;
    double *pcoord;
    int *pconnectivity, *psegment, *psegflag;

    if (NDIMS == 2) {
#ifndef THREED

	/* Define 4 corner points of the rectangle, with this order:
         *            BOUNDZ1
         *          0 ------- 3
         *  BOUNDX0 |         | BOUNDX1
         *          1 ------- 2
         *            BOUNDZ0
         */
        // corner 0
        points[0] = 0;
	points[1] = 0;
        // corner 1
	points[2] = 0;
	points[3] = -param.mesh.zlength;
        // corner 2
	points[4] = param.mesh.xlength;
	points[5] = -param.mesh.zlength;
        // corner 3
	points[6] = param.mesh.xlength;
	points[7] = 0;
	
	for (int i=0; i<n_init_segments; ++i) {
            // 0th node of the i-th segment
	    init_segments[2*i] = i;
            // 1st node of the i-th segment
	    init_segments[2*i+1] = i+1;
	}
	// the 1st node of the last segment is connected back to the 0th node
	init_segments[2*n_init_segments-1] = 0;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDZ0;
        init_segflags[2] = BOUNDX1;
        init_segflags[3] = BOUNDZ1;

        const double max_triangle_size = 1.5 * param.mesh.resolution
            * param.mesh.resolution;

        /********************************************************/
	triangulate_polygon(param.mesh.min_angle, max_triangle_size,
			    npoints, n_init_segments, points,
			    init_segments, init_segflags,
			    &nnode, &nelem, &nseg,
			    &pcoord, &pconnectivity, 
			    &psegment, &psegflag);
        /********************************************************/

        if (nelem <= 0) {
            std::cerr << "Error: triangulation failed\n";
            std::exit(10);
        }
#endif
    } else {
#ifdef THREED
	/* Define 8 corner points of the box, with this order:
         *         4 ------- 7
         *        /         /|
         *       /         / 6
         *      0 ------- 3 /
         *      |         |/
         *      1 ------- 2
         *
         * Cut-out diagram with boundary flag:
         *             4 ------- 7
         *             | BOUNDZ1 |
         *   4 ------- 0 ------- 3 ------- 7 ------- 4
         *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
         *   5 ------- 1 ------- 2 ------- 6 ------- 5
         *             | BOUNDZ0 |
         *             5 ------- 6
         */

        // corner 0
        points[0] = 0;
	points[1] = 0;
	points[2] = 0;
        // corner 1
	points[3] = 0;
	points[4] = 0;
	points[5] = -param.mesh.zlength;
        // corner 2
	points[6] = param.mesh.xlength;
	points[7] = 0;
	points[8] = -param.mesh.zlength;
        // corner 3
	points[9] = param.mesh.xlength;
	points[10] = 0;
	points[11] = 0;
        // corner 4
	points[12] = 0;
	points[13] = param.mesh.ylength;
	points[14] = 0;
        // corner 5
	points[15] = 0;
	points[16] = param.mesh.ylength;
	points[17] = -param.mesh.zlength;
        // corner 6
	points[18] = param.mesh.xlength;
	points[19] = param.mesh.ylength;
	points[20] = -param.mesh.zlength;
        // corner 7
	points[21] = param.mesh.xlength;
	points[22] = param.mesh.ylength;
	points[23] = 0;

        // BOUNDX0
        init_segments[0] = 0;
        init_segments[1] = 1;
        init_segments[2] = 5;
        init_segments[3] = 4;
        // BOUNDY0
        init_segments[4] = 0;
        init_segments[5] = 3;
        init_segments[6] = 2;
        init_segments[7] = 1;
        // BOUNDZ0
        init_segments[8] = 1;
        init_segments[9] = 2;
        init_segments[10] = 6;
        init_segments[11] = 5;
        // BOUNDX1
        init_segments[12] = 3;
        init_segments[13] = 7;
        init_segments[14] = 6;
        init_segments[15] = 2;
        // BOUNDY1
        init_segments[16] = 7;
        init_segments[17] = 4;
        init_segments[18] = 5;
        init_segments[19] = 6;
        // BOUNDZ1
        init_segments[20] = 0;
        init_segments[21] = 4;
        init_segments[22] = 7;
        init_segments[23] = 3;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDY0;
        init_segflags[2] = BOUNDZ0;
        init_segflags[3] = BOUNDX1;
        init_segflags[4] = BOUNDY1;
        init_segflags[5] = BOUNDZ1;

        const double max_tet_size = 0.7 * param.mesh.resolution
            * param.mesh.resolution * param.mesh.resolution;

        /***************************************************************/
	tetrahedralize_polyhedron(param.mesh.max_ratio,
                                  param.mesh.min_tet_angle, max_tet_size,
                                  npoints, n_init_segments, points,
                                  init_segments, init_segflags,
                                  &nnode, &nelem, &nseg,
                                  &pcoord, &pconnectivity,
                                  &psegment, &psegflag);
        /***************************************************************/

#endif
    }

    delete [] points;
    delete [] init_segments;
    delete [] init_segflags;

    var.nnode = nnode;
    var.nelem = nelem;
    var.nseg = nseg;
    var.coord = new arrayd2(pcoord, nnode);
    var.connectivity = new conn_t(pconnectivity, nelem);
    var.segment = new segment_t(psegment, nseg);
    var.segflag = new segflag_t(psegflag, nseg);
}


static void new_mesh_refined_zone(const Param& param, Variables& var)
{
    const Mesh& m = param.mesh;

    // To prevent the meshing library giving us a regular grid, the nodes
    // will be shifted randomly by a small distance
    const double shift_factor = 0.1;

    // typical distance between nodes in the refined zone
    const double d = param.mesh.resolution / std::sqrt(2);

    // adjust the bounds of the refined zone so that nodes are not on the boundary
    double x0, x1, y0, y1, z0, z1;
    double dx, dy, dz;
    int nx, ny, nz;
    x0 = std::max(m.refined_zonex.first, d / m.xlength);
    x1 = std::min(m.refined_zonex.second, 1 - d / m.xlength);
    nx = m.xlength * (x1 - x0) / d + 1;
    dx = m.xlength * (x1 - x0) / nx;
    z0 = std::max(m.refined_zonez.first, d / m.zlength);
    z1 = std::min(m.refined_zonez.second, 1 - d / m.zlength);
    nz = m.zlength * (z1 - z0) / d + 1;
    dz = m.zlength * (z1 - z0) / (nz - 1);

    int npoints;
    if (NDIMS == 2) {
        npoints = nx * nz + 4 * (NDIMS - 1);
    }
    else {
        y0 = std::max(m.refined_zoney.first, d / m.ylength);
        y1 = std::min(m.refined_zoney.second, 1 - d / m.ylength);
        ny = m.ylength * (y1 - y0) / d + 1;
        dy = m.ylength * (y1 - y0) / (ny - 1);
        npoints = nx * ny * nz + 4 * (NDIMS - 1);
    }
    double *points = new double[npoints*NDIMS];

    int n_init_segments = 2 * NDIMS; // 2D:4;  3D:6
    int n_segment_nodes = 2 * (NDIMS - 1); // 2D:2; 3D:4
    int *init_segments = new int[n_init_segments*n_segment_nodes];
    int *init_segflags = new int[n_init_segments];

    int nnode, nelem, nseg;
    double *pcoord;
    int *pconnectivity, *psegment, *psegflag;

    if (NDIMS == 2) {
#ifndef THREED

	/* Define 4 corner points of the rectangle, with this order:
         *            BOUNDZ1
         *          0 ------- 3
         *  BOUNDX0 |         | BOUNDX1
         *          1 ------- 2
         *            BOUNDZ0
         */
        // corner 0
        points[0] = 0;
	points[1] = 0;
        // corner 1
	points[2] = 0;
	points[3] = -param.mesh.zlength;
        // corner 2
	points[4] = param.mesh.xlength;
	points[5] = -param.mesh.zlength;
        // corner 3
	points[6] = param.mesh.xlength;
	points[7] = 0;

        // add refined nodes
        int n = 8;
        for (int i=0; i<nx; ++i) {
            for (int k=0; k<nz; ++k) {
                double rx = drand48() - 0.5;
                double rz = drand48() - 0.5;
                points[n  ] = x0 * m.xlength + (i + shift_factor*rx) * dx;
                points[n+1] = (1-z0) * -m.zlength + (k + shift_factor*rz) * dz;
                n += NDIMS;
            }
        }

	for (int i=0; i<n_init_segments; ++i) {
            // 0th node of the i-th segment
	    init_segments[2*i] = i;
            // 1st node of the i-th segment
	    init_segments[2*i+1] = i+1;
	}
	// the 1st node of the last segment is connected back to the 0th node
	init_segments[2*n_init_segments-1] = 0;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDZ0;
        init_segflags[2] = BOUNDX1;
        init_segflags[3] = BOUNDZ1;

        /********************************************************/
	triangulate_polygon(param.mesh.min_angle, 40*d*d,
			    npoints, n_init_segments, points,
			    init_segments, init_segflags,
			    &nnode, &nelem, &nseg,
			    &pcoord, &pconnectivity,
			    &psegment, &psegflag);
        /********************************************************/

        if (nelem <= 0) {
            std::cerr << "Error: triangulation failed\n";
            std::exit(10);
        }
#endif
    } else {
#ifdef THREED
	/* Define 8 corner points of the box, with this order:
         *         4 ------- 7
         *        /         /|
         *       /         / 6
         *      0 ------- 3 /
         *      |         |/
         *      1 ------- 2
         *
         * Cut-out diagram with boundary flag:
         *             4 ------- 7
         *             | BOUNDZ1 |
         *   4 ------- 0 ------- 3 ------- 7 ------- 4
         *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
         *   5 ------- 1 ------- 2 ------- 6 ------- 5
         *             | BOUNDZ0 |
         *             5 ------- 6
         */

        // corner 0
        points[0] = 0;
	points[1] = 0;
	points[2] = 0;
        // corner 1
	points[3] = 0;
	points[4] = 0;
	points[5] = -param.mesh.zlength;
        // corner 2
	points[6] = param.mesh.xlength;
	points[7] = 0;
	points[8] = -param.mesh.zlength;
        // corner 3
	points[9] = param.mesh.xlength;
	points[10] = 0;
	points[11] = 0;
        // corner 4
	points[12] = 0;
	points[13] = param.mesh.ylength;
	points[14] = 0;
        // corner 5
	points[15] = 0;
	points[16] = param.mesh.ylength;
	points[17] = -param.mesh.zlength;
        // corner 6
	points[18] = param.mesh.xlength;
	points[19] = param.mesh.ylength;
	points[20] = -param.mesh.zlength;
        // corner 7
	points[21] = param.mesh.xlength;
	points[22] = param.mesh.ylength;
	points[23] = 0;

        // add refined nodes
        int n = 24;
        for (int i=0; i<nx; ++i) {
            for (int j=0; j<ny; ++j) {
                for (int k=0; k<nz; ++k) {
                    double rx = drand48() - 0.5;
                    double ry = drand48() - 0.5;
                    double rz = drand48() - 0.5;
                    points[n  ] = x0 * m.xlength + (i + shift_factor*rx) * dx;
                    points[n+1] = y0 * m.ylength + (j + shift_factor*ry) * dy;
                    points[n+2] = (1-z0) * -m.zlength + (k + shift_factor*rz) * dz;
                    n += NDIMS;
                }
            }
        }

        // BOUNDX0
        init_segments[0] = 0;
        init_segments[1] = 1;
        init_segments[2] = 5;
        init_segments[3] = 4;
        // BOUNDY0
        init_segments[4] = 0;
        init_segments[5] = 3;
        init_segments[6] = 2;
        init_segments[7] = 1;
        // BOUNDZ0
        init_segments[8] = 1;
        init_segments[9] = 2;
        init_segments[10] = 6;
        init_segments[11] = 5;
        // BOUNDX1
        init_segments[12] = 3;
        init_segments[13] = 7;
        init_segments[14] = 6;
        init_segments[15] = 2;
        // BOUNDY1
        init_segments[16] = 7;
        init_segments[17] = 4;
        init_segments[18] = 5;
        init_segments[19] = 6;
        // BOUNDZ1
        init_segments[20] = 0;
        init_segments[21] = 4;
        init_segments[22] = 7;
        init_segments[23] = 3;

        // boundary flags (see definition in constants.hpp)
        init_segflags[0] = BOUNDX0;
        init_segflags[1] = BOUNDY0;
        init_segflags[2] = BOUNDZ0;
        init_segflags[3] = BOUNDX1;
        init_segflags[4] = BOUNDY1;
        init_segflags[5] = BOUNDZ1;

        /***************************************************************/
	tetrahedralize_polyhedron(param.mesh.max_ratio,
                                  param.mesh.min_tet_angle, 40*d*d*d,
                                  npoints, n_init_segments, points,
                                  init_segments, init_segflags,
                                  &nnode, &nelem, &nseg,
                                  &pcoord, &pconnectivity,
                                  &psegment, &psegflag);
        /***************************************************************/

#endif
    }

    delete [] points;
    delete [] init_segments;
    delete [] init_segflags;

    var.nnode = nnode;
    var.nelem = nelem;
    var.nseg = nseg;
    var.coord = new arrayd2(pcoord, nnode);
    var.connectivity = new conn_t(pconnectivity, nelem);
    var.segment = new segment_t(psegment, nseg);
    var.segflag = new segflag_t(psegflag, nseg);
}


void create_boundary_flags(Variables& var)
{
    // allocate and init to 0
    var.bcflag = new int_vec(var.nnode);

    // alias for convienence
    int_vec &bcflag = *var.bcflag;
    for (int i=0; i<var.nseg; ++i) {
        int flag = (*var.segflag)[i][0];
        int *n = (*var.segment)[i];
        for (int j=0; j<NDIMS; ++j) {
            bcflag[n[j]] |= flag;
        }
    }
}


void create_boundary_facets(Variables& var)
{
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        for (int i=0; i<FACETS_PER_ELEM; ++i) {
            // set all bits to 1
            int flag = BOUNDX0 | BOUNDX1 | BOUNDY0 | BOUNDY1 | BOUNDZ0 | BOUNDZ1;
            for (int j=0; j<NODES_PER_FACET; ++j) {
                // find common flags
                int n = NODE_OF_FACET[i][j];
                flag &= (*var.bcflag)[conn[n]];
            }
            if (flag) {
                // this facet belongs to a boundary
                int n = bdry_order.find(flag)->second;
                var.bfacets[n].push_back(std::make_pair(e, i));
            }
        }
    }

    // for (int n=0; n<6; ++n) {
    //     std::cout << "boundary facet " << n << ":\n";
    //     print(std::cout, var.bfacets[n]);
    //     std::cout << '\n';
    //     for (int i=0; i<var.bfacets[n].size(); ++i) {
    //         int e = var.bfacets[n][i].first;
    //         int f = var.bfacets[n][i].second;
    //         const int *conn = &(*var.connectivity)[e][0];
    //         std::cout << i << ", " << e << ":";
    //         for (int j=0; j<NODES_PER_FACET; ++j) {
    //             std::cout << " " << conn[NODE_OF_FACET[f][j]];
    //         }
    //         std::cout << '\n';
    //     }
    //     std::cout << '\n';
    // }
}


void create_support(Variables& var)
{
    var.support = new std::vector<int_vec>(var.nnode);

    // create the inverse mapping of connectivity
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            (*var.support)[conn[i]].push_back(e);
        }
    }
    // std::cout << "support:\n";
    // print(std::cout, *var.support);
    // std::cout << "\n";
}


void create_elem_groups(Variables& var)
{
    /* Decompose the whole mesh into groups of elements, each element
     * is the member of one and only one group. The elements in one group
     * are disjoint, i.e. not sharing the same nodes, edges nor faces.
     *
     * These groups are the boundary of openmp parallelism. Within each
     * group, it is safe to use openmp.
     */

    var.egroups = new std::vector<int_vec>;

#ifdef USE_OMP
    // element # ordered by "crowdness" (see below)
    std::vector<std::size_t> crowdest_elem(var.nelem);
    {
        // how many elements are sharing nodes with this element?
        std::vector<std::size_t> crowdness(var.nelem);
        for (int e=0; e<var.nelem; ++e) {
            const int *conn = (*var.connectivity)[e];
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                int n = conn[i];
                crowdness[e] += (*var.support)[n].size();
            }
        }

        sortindex_reversed(crowdness, crowdest_elem);
        // print(std::cout, crowdness);
        // std::cout << "\n\n";
        // print(std::cout, crowdest_elem);
        // std::cout << "\n\n";
        // for (int e=0; e<var.nelem; ++e) {
        //     std::cout << e << " : " << crowdness[crowdest_elem[e]] << ", " << crowdest_elem[e] << '\n';
        // }
    }

    int start = 0;
    const std::size_t sentinel = crowdest_elem[0];
    while (1) {
        // gp is the current egroup
        var.egroups->push_back(int_vec());
        int_vec& gp = var.egroups->back();

        // nodes taken by the current egroup so far
        std::vector<bool> nodes_taken(var.nnode, 0);

        // the starting element is always available to take
        {
            int e = crowdest_elem[start];
            const int *conn = (*var.connectivity)[e];
            // mark nodes as taken
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                int n = conn[i];
                nodes_taken[n] = 1;
            }
            gp.push_back(e);
            // mark element as taken
            crowdest_elem[start] = sentinel;
        }
        for (int ee=start+1; ee<var.nelem; ++ee) {
            std::size_t e = crowdest_elem[ee];
            // this element is taken, skip
            if (e == sentinel) continue;

            const int *conn = (*var.connectivity)[e];

            // does this element share any node with other elements in the group?
            bool is_sharing_node = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                int n = conn[i];
                if (nodes_taken[n]) {
                    is_sharing_node = 1;
                    break;
                }
            }
            // skip this element
            if (is_sharing_node) continue;

            // None of the nodes are taken. This element belongs to this group
            // mark nodes as taken
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                int n = conn[i];
                nodes_taken[n] = 1;
            }
            gp.push_back(e);
            // mark element as taken
            crowdest_elem[ee] = sentinel;
        }

        // before starting over, find the first available element
        start++;
        bool found = 0;
        for (int ee=start; ee<var.nelem; ++ee) {
            if (crowdest_elem[ee] != sentinel) {
                found = 1;
                start = ee;
                break;
            }
        }

        // none of the elements are available, job is done
        if (! found) break;
    }

#else

    // Not using openmp, only need one group for all elements
    var.egroups->push_back(int_vec(var.nelem));
    int_vec& gp = var.egroups->back();
    for (int e=0; e<var.nelem; ++e) {
        gp[e] = e;
    }

#endif

    // std::cout << "egroups:" << var.egroups->size() << " groups\n";
    // for (int i=0; i<var.egroups->size(); i++) {
    //     std::cout << (*var.egroups)[i].size() << '\n';
    // }

    // print(std::cout, *var.egroups);
    // std::cout << '\n';
}


void create_new_mesh(const Param& param, Variables& var)
{
    switch (param.mesh.meshing_option) {
    case 1:
        new_mesh_uniform_resolution(param, var);
        break;
    case 2:
        new_mesh_refined_zone(param, var);
        break;
    default:
        std::cout << "Error: unknown meshing option: " << param.mesh.meshing_option << '\n';
        std::exit(1);
    }
    // std::cout << "segment:\n";
    // print(std::cout, *var.segment);
    // std::cout << '\n';
    // std::cout << "segflag:\n";
    // print(std::cout, *var.segflag);
    // std::cout << '\n';

    create_boundary_flags(var);
    create_boundary_facets(var);
    create_support(var);
    create_elem_groups(var);
}
