#ifndef DYNEARTHSOL3D_MARKERSET_HPP
#define DYNEARTHSOL3D_MARKERSET_HPP

class MarkerSet
{

public:
    MarkerSet( const Param& param, Variables& var );
    ~MarkerSet() 
    { 
        delete _id;
        delete _eta; 
        delete _elem;
        delete _mattype;
    }

    void write();
    void read();
    void resize(const int);

    inline int get_nmarkers() { return _nmarkers; }
    inline void set_nmarkers(int n) { _nmarkers = n; }

    inline int get_id(int m) { return (*_id)[m]; }
    inline void set_id(const int m, const int i) { (*_id)[m] = i; }
    inline void resize_id(const int n) { (*_id).resize(n); }

    inline int get_elem(int m) { return (*_elem)[m]; }
    inline void set_elem(const int m, const int e) { (*_elem)[m] = e; }
    inline void resize_elem(const int n) { (*_elem).resize(n); }

    inline int get_mattype(int m) { return (*_mattype)[m]; }
    inline void set_mattype(const int m, const int mt) { (*_mattype)[m] = mt; }
    inline void resize_mattype(const int n) { (*_mattype).resize(n); }

    inline double *get_eta(int m) { return (*_eta)[m]; }
    inline void set_eta( const int i, const double r[NDIMS] );
    inline void resize_eta( double *new_eta, const int n ) { (*_eta).reset( new_eta, n); }

private:
    // Didn't create a data type for an individual marker to follow the "structure of arrays" concept.
    // Number of markers (may change during remeshing)
    int _nmarkers;
    // Barycentric (local) coordinate within the reference element
    shapefn *_eta;
    // Containing element
    int_vec *_elem;
    // Material type
    int_vec *_mattype;
    // Unique id
    int_vec *_id;

    void random_markers( const Param&, Variables& );
    void allocate_markerdata( const int, const int );

};

void remap_markers(const Param&, Variables &, 
                   const array_t &, const conn_t &);

#endif
