#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
#include "binaryio.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "geometry.hpp"
#include "utils.hpp"

namespace {

    const int DEBUG = 0;
    const double over_alloc_ratio = 2.0;  // how many extra space to allocate for future expansion

}


MarkerSet::MarkerSet(const Param& param, Variables& var)
{
    _last_id = _nmarkers = 0;

    switch ( param.markers.init_marker_option ) {
    case 1:
        random_markers(param, var);
        break;
    default:
        std::cerr << "Error: unknown init_marker_option: " << param.markers.init_marker_option << ". The only valid option is '1'.\n";
        break;
    }
}


MarkerSet::MarkerSet(const Param& param, Variables& var, BinaryInput& bin)
{
    // init from checkpoint file
    read(var, bin);
}


void MarkerSet::allocate_markerdata( const int max_markers )
{
    _reserved_space = max_markers;
    _eta = new shapefn( max_markers );
    _elem = new int_vec( max_markers );
    _mattype = new int_vec( max_markers );
    _id = new int_vec( max_markers );
}


void MarkerSet::random_eta( double *eta )
{
    // eta for randomly scattered markers within an element.
    // An alternative would be to fix barycentric coordinates and add random perturbations.
    //
    // 1. populate eta with random numbers between 0 and 1.0.
    double sum = 0;
    for( int n = 0; n < NODES_PER_ELEM; n++ ) {
        eta[n] = (rand()/(double)RAND_MAX);
        sum += eta[n];
    }
    // 2. normalize.
    double inv_sum = 1.0/sum;
    for( int n = 0; n < NODES_PER_ELEM; n++ )
        eta[n] *= inv_sum;
}


void MarkerSet::append_marker( double *eta, int el, int mt )
{
    // Ensure sufficient array size
    if( _nmarkers == _reserved_space ) {
        // Resize the marker-related arrays if necessary.
        const int newsize = _nmarkers * over_alloc_ratio;
        resize( newsize );
    }

    int m = _nmarkers;
    std::memcpy((*_eta)[m], eta, NODES_PER_ELEM*sizeof(double));
    (*_elem)[m] = el;
    (*_mattype)[m] = mt;
    (*_id)[m] = _last_id;

    if(DEBUG > 1) {
        std::cout << el << " " << m << " "
                  << _nmarkers << " " << (*_mattype)[_nmarkers] << " "
                  << eta[0] << "+" << eta[1] << "+" << eta[2];
#ifdef THREED
        std::cout << "+" << eta[3]
                  << "=" << (eta[0]+eta[1]+eta[2]+eta[3]) << "\n";
#else
        std::cout << "=" << (eta[0]+eta[1]+eta[2]) << "\n";
#endif
    }

    ++_nmarkers;
    ++_last_id;
}


void MarkerSet::append_random_marker_in_elem( int el, int mt )
{
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt);
}


void MarkerSet::random_markers( const Param& param, Variables &var )
{
    const int ne = var.nelem;
    const int mpe = param.markers.markers_per_element;
    const int num_markers = ne*mpe;
    const int max_markers = num_markers * over_alloc_ratio;

    // allocate memory for data members.
    allocate_markerdata( max_markers );
    
    // initialize random seed:
    srand (time(NULL));
    
    // Generate particles in each element.
    for( int e = 0; e < ne; e++ )
        for( int m = 0; m < mpe; m++ ) {
            // random barycentric coordinate
            double eta[NODES_PER_ELEM];
            random_eta(eta);

            // decide the mattype of markers
            int mt;
            if (param.ic.mattype_option == 0) {
                mt = (int)(*((*var.regattr)[e])); // mattype should take a reginal attribute assigned during meshing.
            }
            else {
                double p[NDIMS] = {0};
                const int *conn = (*var.connectivity)[e];
                for(int i=0; i<NDIMS; i++) {
                    for(int j=0; j<NODES_PER_ELEM; j++)
                        p[i] += (*var.coord)[ conn[j] ][i] * eta[j];
                }
                // modify mt according to the marker coordinate p
                switch (param.ic.mattype_option) {
                case 1:
                    // lower half: 1; upper half: 0
                    if (p[NDIMS-1] < -0.5 * param.mesh.zlength)
                        mt = 1;
                    else
                        mt = 0;
                    break;
                default:
                    std::cerr << "Error: unknown ic.mattype_option: " << param.ic.mattype_option << '\n';
                    std::exit(1);
                }
            }
            append_marker(eta, e, mt);
            ++(*var.elemmarkers)[e][mt];
        }
}


void MarkerSet::remove_marker(int i)
{
    // Replace marker i by the last marker.
    --_nmarkers;
    std::memcpy( (*_eta)[i], (*_eta)[_nmarkers], sizeof(double)*(NODES_PER_ELEM) );
    (*_id)[i] = (*_id)[_nmarkers];
    (*_elem)[i] = (*_elem)[_nmarkers];
    (*_mattype)[i] = (*_mattype)[_nmarkers];
}

    
void MarkerSet::set_eta( const int i, const double r[NDIMS] ) {
    double sum = 0.0;
    for( int j = 0; j < NDIMS; j++ ) {
        (*_eta)[i][j] = r[j];
        sum += r[j];
    }
    (*_eta)[i][NDIMS] = 1.0 - sum;
}


void MarkerSet::resize( const int newsize )
{
    if( newsize > _reserved_space ) {
        // enlarge arrays
        std::cout << "  Increasing marker arrays size to " << newsize << " markers.\n";
        _reserved_space = newsize;

        double *new_eta = new double[NODES_PER_ELEM * newsize];
        std::copy( (*_eta)[0], (*_eta)[_nmarkers], new_eta );
        _eta->reset( new_eta, _nmarkers );

        _elem->resize( newsize );
        _mattype->resize( newsize );
        _id->resize( newsize );
    }
    // else if( nmarkers_new < _reserved_space ) {
    //     // TBD: shrink arrays
    //     _reserved_space = newsize;
    // }
    else {
        // New size is too close to old size, don't do anything.
    }
}


void remap_markers(const Param& param, Variables &var, const array_t &old_coord, 
                   const conn_t &old_connectivity)
{
    // Re-create elemmarkers
    delete var.elemmarkers;
    create_elemmarkers( param, var );

    double_vec new_volume( var.nelem );
    compute_volume( *var.coord, *var.connectivity, new_volume );

    Barycentric_transformation bary( *var.coord, *var.connectivity, new_volume );

    // nearest-neighbor search structure
    double **centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements
    ANNkd_tree kdtree(centroid, var.nelem, NDIMS);
    const int k = std::min(20, var.nelem);  // how many nearest neighbors to search?
    const double eps = 0.001 * param.mesh.resolution; // tolerance of distance error
    int *nn_idx = new int[k];
    double *dd = new double[k];

    // Loop over all the old markers and identify a containing element in the new mesh.
    MarkerSet *ms = var.markerset; // alias to var.markerset
    int last_marker = ms->get_nmarkers();
    int i = 0;
    while (i < last_marker) {
        bool found = false;

        // 1. Get physical coordinates, x, of an old marker.
        int eold = ms->get_elem(i);
        double x[NDIMS] = {0};
        for (int j = 0; j < NDIMS; j++)
            for (int k = 0; k < NODES_PER_ELEM; k++)
                x[j] += ms->get_eta(i)[k]*
                    old_coord[ old_connectivity[eold][k] ][j];

        if (DEBUG) {
            std::cout << "marker #" << i << " old_elem " << eold << " x: ";
            print(std::cout, x, NDIMS);
        }

        // 2. look for nearby elements.
        // Note: kdtree.annkSearch() is not thread-safe, cannot use openmp in this loop
        kdtree.annkSearch(x, k, nn_idx, dd, eps);

        for( int j = 0; j < k; j++ ) {
            int e = nn_idx[j];
            double r[NDIMS];

            bary.transform(x, e, r);

            // change this if-condition to (i == N) to debug the N-th marker
            if (0) {
                std::cout << '\n' << j << " check elem #" << e << ' ';
                print(std::cout, r, NDIMS);
            }

            if (bary.is_inside(r)) {
                ms->set_eta(i, r);
                ms->set_elem(i, e);
                ++(*(var.elemmarkers))[e][ms->get_mattype(i)];
            
                found = true;
                ++i;
                if (DEBUG) {
                    std::cout << " in element " << e << '\n';
                }
                break;
            }
        }

        if( found ) continue;

        if (DEBUG) {
            std::cout << " not in any element" << '\n';
        }

        /* not found */
        {
            // Since no containing element has been found, delete this marker.
            // Note i is not inc'd.
            --last_marker;
            ms->remove_marker(i);
        }
    }

    delete [] nn_idx;
    delete [] dd;
    delete [] centroid[0];
    delete [] centroid;

    // If any new element has too few markers, generate markers in them.
    const int mpe = param.markers.markers_per_element;
    for( int e = 0; e < var.nelem; e++ ) {
        int num_marker_in_elem = 0;
        for( int i = 0; i < param.mat.nmat; i++ )
            num_marker_in_elem += (*(var.elemmarkers))[e][i];

        if (num_marker_in_elem < mpe / 2) {  // mpe must >= 2
            // cumulative probability density function of mattype in this element
            double_vec cpdf(param.mat.nmat);
            cpdf[0] = (*(var.elemmarkers))[e][0] / double(num_marker_in_elem);
            for( int i = 1; i < param.mat.nmat - 1; i++ )
                cpdf[i] = cpdf[i-1] + (*(var.elemmarkers))[e][i] / double(num_marker_in_elem);
            cpdf[param.mat.nmat - 1] = 1; // fix to 1 to avoid round-off error
            while( num_marker_in_elem < mpe / 2 ) {
                // Determine new marker's matttype based on cpdf
                auto upper = std::upper_bound(cpdf.begin(), cpdf.end(), drand48());
                const int mt = upper - cpdf.begin();
                ms->append_random_marker_in_elem(e, mt);

                ++(*var.elemmarkers)[e][mt];
                ++num_marker_in_elem;
            }
        }
    }
}


void MarkerSet::write(BinaryOutput &bin)
{
    int_vec itmp(2);
    itmp[0] = _nmarkers;
    itmp[1] = _last_id;
    bin.write_array(itmp, "markerset size");

    bin.write_array(*_eta, "markerset.eta");
    bin.write_array(*_elem, "markerset.elem");
    bin.write_array(*_mattype, "markerset.mattype");
    bin.write_array(*_id, "markerset.id");

}


void MarkerSet::read(Variables &var, BinaryInput &bin)
{
    int_vec itmp(2);
    bin.read_array(itmp, "markerset size");
    _nmarkers = itmp[0];
    _last_id = itmp[1];

    allocate_markerdata(_nmarkers);

    bin.read_array(*_eta, "markerset.eta");
    bin.read_array(*_elem, "markerset.elem");
    bin.read_array(*_mattype, "markerset.mattype");
    bin.read_array(*_id, "markerset.id");

    for( int i = 0; i < _nmarkers; i++ ) {
        int e = (*_elem)[i];
        int mt = (*_mattype)[i];
        ++(*var.elemmarkers)[e][mt];
    }
}
