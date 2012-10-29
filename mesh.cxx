
#include <cstdio>
#include <cstring>

#ifdef THREED

#include "tetgen/tetgen.h"

#else

#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS

#endif //ifdef THREED

#include "constants.hpp"
#include "parameters.hpp"
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
    struct triangulateio in, out;

    // add 'Q' for no output; add multiple 'V's for verbose output
    std::sprintf(options, "pq%fjza%f", min_angle, max_area);
    std::puts(options);

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
(double min_angle, double max_volume,
 int npoints, int nsegments,
 double *points, int *segments, int *segflags,
 int *noutpoints, int *ntriangles, int *noutsegments,
 double **outpoints, int **triangles,
 int **outsegments, int **outsegflags,
 tetgenio **out)
{
#ifdef THREED
    //
    // Setting Tetgen options.
    // Using "tetgenhavior", rather than string, to specify the options for its flexibility.
    // However, we have to maintain consistency between these options.
    // See tetgenbehavior::parse_commandline() for details.
    //
    tetgenbehavior options;

    options.quiet = 0;
    options.verbose = 1;

    options.plc = 1;
    options.zeroindex = 0;
    options.vtkview = 1;

    options.quality = 1;
    options.minratio = 1.44;  // TODO: provide an input parameter for tuning it
    options.mindihedral = min_angle;
    options.maxdihedral = 180 - 3 * options.mindihedral;
    options.fixedvolume = 1;
    options.maxvolume = max_volume;

    // TODO: to be replaced by param.modelname
    std::strncpy(options.outfilename, "zzz", tetgenio::FILENAMESIZE-1);

    // derived from options above
    options.useshelles = 1;
    options.goodratio = options.minratio;
    options.goodratio *= options.goodratio;

    //
    // Setting input arrays to tetgen
    //
    tetgenio *in;
    in = new tetgenio;
    in->pointlist = points;
    in->numberofpoints = npoints;

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

    in->facetlist = fl;
    in->facetmarkerlist = segflags;
    in->numberoffacets = nsegments;

    in->holelist = NULL;
    in->numberofholes = 0;

    in->regionlist = NULL;
    in->numberofregions = 0;

    *out = new tetgenio;
    /*******************************/
    tetrahedralize(&options, in, *out);
    /*******************************/

    // the destructor of tetgenio will free any non-NULL pointer
    // set in.pointers to NULL to prevent double-free
    in->pointlist = NULL;
    in->facetmarkerlist = NULL;
    in->facetlist = NULL;
    delete [] polys;
    delete [] fl;
    delete in;

    *noutpoints = (*out)->numberofpoints;
    *outpoints = (*out)->pointlist;

    *ntriangles = (*out)->numberoftetrahedra;
    *triangles = (*out)->tetrahedronlist;

    *noutsegments = (*out)->numberoftrifaces;
    *outsegments = (*out)->trifacelist;
    *outsegflags = (*out)->trifacemarkerlist;

#endif
}


void create_boundary(const Param& param, Variables& var)
{
    // allocate and init to 0
    var.bcflag = new int_vec(var.nnode);

    // alias for convienence
    int_vec &bcflag = *var.bcflag;
    for (int i=0; i<var.nseg; ++i) {
        int flag = (*var.segflag)[i];
        int *n = &(*var.segment)[i][0];
        bcflag[n[0]] |= flag;
        bcflag[n[1]] |= flag;
    }
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

        const double min_triangle_angle = 32.;
        const double max_triangle_size = 1.5 * param.mesh.resolution
            * param.mesh.resolution;

	triangulate_polygon(min_triangle_angle, max_triangle_size,
			    npoints, n_init_segments, points,
			    init_segments, init_segflags,
			    &nnode, &nelem, &nseg,
			    &pcoord, &pconnectivity, 
			    &psegment, &psegflag);

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

        const double min_tet_angle = 22.;
        const double max_tet_size = 0.7 * param.mesh.resolution
            * param.mesh.resolution * param.mesh.resolution;
	tetrahedralize_polyhedron(min_tet_angle, max_tet_size,
                                  npoints, n_init_segments, points,
                                  init_segments, init_segflags,
                                  &nnode, &nelem, &nseg,
                                  &pcoord, &pconnectivity,
                                  &psegment, &psegflag, &var.tetgen);
#endif
    }

    delete [] points;
    delete [] init_segments;
    delete [] init_segflags;

    var.nnode = nnode;
    var.nelem = nelem;
    var.nseg = nseg;
    var.coord = new double2d_ref(pcoord, boost::extents[nnode][NDIMS]);
    var.connectivity = new int2d_ref(pconnectivity, boost::extents[nelem][NODES_PER_ELEM]);
    var.segment = new int2d_ref(psegment, boost::extents[nseg][NDIMS]);
    var.segflag = new int1d_ref(psegflag, boost::extents[nseg]);

    /*
    std::cout << "coord:\n";
    print(std::cout, *var.coord);
    std::cout << "\n";

    std::cout << "connectivity:\n";
    print(std::cout, *var.connectivity);
    std::cout << "\n";

    std::cout << "segment:\n";
    print(std::cout, *var.segment);
    std::cout << "\n";

    std::cout << "segflag:\n";
    print(std::cout, *var.segflag);
    std::cout << "\n";
    */
}



void create_new_mesh(const Param& param, Variables& var)
{
    new_mesh_uniform_resolution(param, var);
    create_boundary(param, var);
}
