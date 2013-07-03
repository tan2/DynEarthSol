#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
#include "markerset.hpp"
#include "mesh.hpp"

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
    _nmarkers = max_markers;
    _eta = new shapefn( max_markers );
    _elem = new int_vec( max_markers );
    _mattype = new int_vec( max_markers );
    _id = new int_vec( max_markers );
}


void MarkerSet::random_markers( const Param& param, Variables &var )
{
    const int ne = var.nelem;
    const int mpe = param.markers.markers_per_element;
    const int max_markers = ne*mpe;
    
    // allocate memory for data members.
    allocate_markerdata( max_markers, mpe );
    
    // initialize random seed:
    srand (time(NULL));
    
    // Store the number of markers
    _nmarkers = max_markers;

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

            if( NDIMS==3 )
                std::cerr << e <<"/"<< ne <<" "<< m <<" "
                          <<pid<<" "<<(*_mattype)[pid]<<" "
                          <<inv_sum<<" "
                          <<" size of eta= "<<(*_eta).size()<<" "
                          <<eta[0]<<"+"<<eta[1]<<"+"<<eta[2]
                          <<"+"<<eta[3]<<"="<<(eta[0]+eta[1]+eta[2]+eta[3])<<"\n";
            else
                std::cerr << e <<"/"<< ne <<" "<< m <<" "
                          <<pid<<" "<<(*_mattype)[pid]<<" "
                          <<inv_sum<<" "
                          <<" size of eta= "<<(*_eta).size()<<" "
                          <<eta[0]<<"+"<<eta[1]<<"+"<<eta[2]
                          <<"="<<(eta[0]+eta[1]+eta[2])<<"\n";
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

void MarkerSet::resize( const int nmarkers_new ) {
    
    double *new_eta = new double[NODES_PER_ELEM * nmarkers_new];
    std::copy( get_eta(0), get_eta(nmarkers_new), new_eta );
    resize_eta( new_eta, nmarkers_new );
    
    resize_id( nmarkers_new );
    resize_elem( nmarkers_new );
    resize_mattype( nmarkers_new );
}

void remap_markers(const Param& param, Variables &var, const array_t &old_coord, 
                   const conn_t &old_connectivity,
                   const Barycentric_transformation &bary)
{

    // Re-create elemmarkers
    delete var.elemmarkers;
    create_elemmarkers( param, var );

    // Loop over all the old markers and identify a containing element in the new mesh.
    MarkerSet *ms = var.markerset; // alias to var.markerset
    const int nmarkers_old = ms->get_nmarkers();
    for(int i = 0 ; i < nmarkers_old ; i++) {

        bool found = false;

        // 1. Get physical coordinates, x, of an old marker.
        double x[NDIMS];
        std::memset( x, 0.0, sizeof(double)*NDIMS );
        for (int j = 0; j < NDIMS; j++)
            for (int k = 0; k < NDIMS+1; k++)
                x[j] += ms->get_eta(i)[k]*
                    old_coord[old_connectivity[ms->get_elem(i)][k]][j];

        // 2. First look into a new element with id = _elem[i].
        {
            double r[NDIMS];
            int e = ms->get_elem(i);
            
            bary.transform(x, e, r);
            if (bary.is_inside(r)) {
                ms->set_eta( i, r );
                ++(*(var.elemmarkers))[e][ms->get_mattype(i)];
                
                found = true;
            }
        }
        if( found ) continue;

        // TODO: usiing kd-tree to find the element
        // 3. If not found, loop over the entire element set.
        for( int e = 0; e < var.nelem && e != ms->get_elem(i); e++ ) {
            double r[NDIMS];

            bary.transform(x, e, r);
            if (bary.is_inside(r)) {
                ms->set_eta(i, r);
                ms->set_elem(i, e);
                ++(*(var.elemmarkers))[e][ms->get_mattype(i)];
            
                found = true;
                break;
            }
        }
        if( found ) continue;
        
        // Since no containing element has been found, delete this marker.
        if( !found ) {
            int last_marker = nmarkers_old-1;

            std::memcpy( ms->get_eta(i), ms->get_eta(last_marker), sizeof(double)*(NODES_PER_ELEM) );
            ms->set_id( i, ms->get_id(last_marker) );
            ms->set_elem( i, ms->get_elem(last_marker) );
            ms->set_mattype( i, ms->get_mattype(last_marker) );

            ms->set_nmarkers( last_marker );
        }
    }
    // Resize the marker-related arrays.
    const int nmarkers_new = ms->get_nmarkers();
    ms->resize( nmarkers_new );

    // TBD: If any new element has too few markers, generate markers in them. 
    const int mpe = param.markers.markers_per_element;
    for( int e = 0; e < var.nelem; e++ ) {
        int num_marker = 0;
        for( int i = 0; i < param.mat.nmat; i++ )
            num_marker += (*(var.elemmarkers))[e][i];
        
        while( num_marker < mpe ) {
            // add a marker: NN sufficient for mattype decision?
            num_marker++;
        }
    }

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
    
}
