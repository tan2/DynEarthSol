
#include <cstdio>
#include <cstdlib>

#include "tetgen/tetgen.h"

#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle/triangle.h"
#undef REAL
#undef VOID
#undef ANSI_DECLARATORS

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
    return;
}


void create_boundary(const Param& param, Variables& var)
{
    // allocate and init to 0
    var.bcflag = new int2d(boost::extents[var.nnode][NDIMS]);
    std::fill(var.bcflag->origin(), var.bcflag->origin() + var.bcflag->size(), 0);

    int2d& bcflag = *var.bcflag;
    for (int i=0; i<var.nseg; ++i) {
        int flag = (*var.segflag)[i];
        switch (flag) {
        case BOUNDX0:
        case BOUNDX1:
            bcflag[i][0] = flag;
            break;
        case BOUNDY0:
        case BOUNDY1:
            if (NDIMS == 2) {
                std::cerr << "Error: this segment flag only works for 3D!\n";
                std::exit(1);
            }
            bcflag[i][1] = flag;
            break;
        case BOUNDZ0:
        case BOUNDZ1:
            bcflag[i][NDIMS-1] = flag;
            break;
        default:
            std::cerr << "Error: Unknown segment flag!\n";
            std::exit(1);
        }
    }
}


static void new_mesh_uniform_resolution(const Param& param, Variables& var)
{
    const double min_triangle_angle = 32.;
    const double max_triangle_size = 1.5 * param.mesh.resolution 
        * param.mesh.resolution;

    int n_init_segments = 4;
    int npoints = 4;
    double *points = (double*) std::malloc(npoints*NDIMS*sizeof(double));
    int *init_segments = (int*) std::malloc(n_init_segments*NDIMS
					    *sizeof(int));
    int *init_segflags = (int*) std::malloc(n_init_segments
					    *sizeof(int));
 
    int nnode, nelem, nseg;
    double *pcoord;
    int *pconnectivity, *psegment, *psegflag;

    if (NDIMS == 2) {

	/* Define 4 corner points of the box, with this order:
         *            BOUNDZ1
         *          0 ------- 3
         *  BOUNDX0 |         | BOUNDX1
         *          1 ------- 2
         *            BOUNDZ0
         */
        points[0] = 0;
	points[1] = 0;
	points[2] = 0;
	points[3] = -param.mesh.zlength;
	points[4] = param.mesh.xlength;
	points[5] = -param.mesh.zlength;
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

        free(points);
	free(init_segments);
	free(init_segflags);

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

    } else {

    }

}



void create_new_mesh(const Param& param, Variables& var)
{
    new_mesh_uniform_resolution(param, var);
    create_boundary(param, var);
}
