#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
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
    switch ( param.markers.init_marker_option ) {
    case 1:
        random_markers(param, var);
        break;
    default:
        std::cerr << "Error: unknown init_marker_option: " << param.markers.init_marker_option << ". The only valid option is '1'.\n";
        break;
    }
}


void MarkerSet::allocate_markerdata( const int max_markers, const int mpe )
{
    _reserved_space = max_markers;
    _eta = new shapefn( max_markers );
    _elem = new int_vec( max_markers );
    _mattype = new int_vec( max_markers );
    _id = new int_vec( max_markers );
}


void MarkerSet::random_markers( const Param& param, Variables &var )
{
    const int ne = var.nelem;
    const int mpe = param.markers.markers_per_element;
    const int num_markers = ne*mpe;
    const int max_markers = num_markers * over_alloc_ratio;

    // allocate memory for data members.
    allocate_markerdata( max_markers, mpe );
    
    // initialize random seed:
    srand (time(NULL));
    
    // Store the number of markers
    _nmarkers = num_markers;

    // Generate particles in each element.
    for( int e = 0; e < ne; e++ )
        for( int m = 0; m < mpe; m++ ) {
            int pid = m + e*mpe;
            
            (*_id)[pid] = pid; 
            (*_elem)[pid] = e;
            (*_mattype)[pid] = (int)(*((*var.regattr)[e])); // mattype should take a reginal attribute assigned during meshing.

            ++(*var.elemmarkers)[e][(*_mattype)[pid]];


            // _eta for randomly scattered markers within an element. 
            // An alternative would be to fix barycentric coordinates and add random perturbations.
            //
            // 1. populate _eta[pid] with random numbers between 0 and 1.0.
            double sum = 0.0;
            double *eta = (*_eta)[pid];
            for( int n = 0; n < NODES_PER_ELEM; n++ ) {
                eta[n] = (rand()/(double)RAND_MAX);
                sum += eta[n];
                // std::cerr << e <<" "<< m <<" "<<pid<<" "<<n<<" "
                //           <<(*_eta)[pid][n]<<" "<<sum<<"\n";
            }
            assert(sum > 0.0);

            // 2. normalize.
            double inv_sum = 1.0/sum;
            for( int n = 0; n < NODES_PER_ELEM; n++ )
                eta[n] *= inv_sum;

            if(DEBUG > 1) {
                std::cout << e << "/" << ne << " " << m << " "
                          << pid << " " << (*_mattype)[pid] << " "
                          << inv_sum << " "
                          << " size(eta) = " << (*_eta).size() << " "
                          << eta[0] << "+" << eta[1] << "+" << eta[2];
#ifdef THREED
                std::cout << "+" << eta[3]
                          << "=" << (eta[0]+eta[1]+eta[2]+eta[3]) << "\n";
#else
                std::cout << "=" << (eta[0]+eta[1]+eta[2]) << "\n";
#endif
            }
        }
}


    
void MarkerSet::write()
{

}

void MarkerSet::read()
{

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
            // Since no containing element has been found, delete this marker,
            // replace it by the last marker. Note i is not inc'd.
            --last_marker;

            std::memcpy( ms->get_eta(i), ms->get_eta(last_marker), sizeof(double)*(NODES_PER_ELEM) );
            ms->set_id( i, ms->get_id(last_marker) );
            ms->set_elem( i, ms->get_elem(last_marker) );
            ms->set_mattype( i, ms->get_mattype(last_marker) );

            ms->set_nmarkers( last_marker );
        }
    }

    delete [] nn_idx;
    delete [] dd;
    delete [] centroid[0];
    delete [] centroid;

    // Resize the marker-related arrays if necessary.
    const int nmarkers_new = ms->get_nmarkers() * over_alloc_ratio;
    ms->resize( nmarkers_new );

    // TBD: If any new element has too few markers, generate markers in them. 
    const int mpe = param.markers.markers_per_element;
    for( int e = 0; e < var.nelem; e++ ) {
        int num_marker = 0;
        for( int i = 0; i < param.mat.nmat; i++ )
            num_marker += (*(var.elemmarkers))[e][i];

        // temporary safeguard, remove it when TBD is finished.
        if( num_marker == 0 ) {
            std::cerr << "Error: no marker in element #" << e << '\n';
            print(std::cout, *var.elemmarkers);
            std::cout << '\n';
            std::exit(10);
        }

        while( num_marker < mpe ) {
            // add a marker: NN sufficient for mattype decision?
            num_marker++;
        }
    }
    /*
    // debug print
    if( NDIMS == 3 )
        for( int i = 0; i < ms->get_nmarkers(); i++ )
            std::cerr << i <<" / "<< ms->get_nmarkers() <<" id="<<ms->get_mattype(i)
                      <<" el="<<ms->get_elem(i)<<" mat="<<ms->get_mattype(i)<<" eta="
                      <<ms->get_eta(i)[0]<<"+"<<ms->get_eta(i)[1]<<"+"
                      <<ms->get_eta(i)[2]<<"+"<<ms->get_eta(i)[3]<<"="
                      <<(ms->get_eta(i)[0]+ms->get_eta(i)[1]+ms->get_eta(i)[2]+ms->get_eta(i)[3])<<"\n";
    else
        for( int i = 0; i < ms->get_nmarkers(); i++ )
            std::cerr << i <<" / "<< ms->get_nmarkers() <<" id="<<ms->get_mattype(i)
                      <<" el="<<ms->get_elem(i)<<" mat="<<ms->get_mattype(i)<<" eta="
                      <<ms->get_eta(i)[0]<<"+"<<ms->get_eta(i)[1]<<"+"
                      <<ms->get_eta(i)[2]<<"="
                      <<(ms->get_eta(i)[0]+ms->get_eta(i)[1]+ms->get_eta(i)[2])<<"\n";
    */
}
