#ifndef DYNEARTHSOL3D_MARKERSET_HPP
#define DYNEARTHSOL3D_MARKERSET_HPP

// forward declaration
class BinaryOutput;

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

    void append_random_marker_in_elem( int el, int mt );
    void remove_marker(int i);
    void resize(const int);
    void write(BinaryOutput &bin);
    void read(BinaryOutput &bin);

    inline int get_nmarkers() const { return _nmarkers; }
    inline void set_nmarkers(int n) { _nmarkers = n; }

    inline int get_id(int m) const { return (*_id)[m]; }
    inline void set_id(const int m, const int i) { (*_id)[m] = i; }

    inline int get_elem(int m) const { return (*_elem)[m]; }
    inline void set_elem(const int m, const int e) { (*_elem)[m] = e; }

    inline int get_mattype(int m) const { return (*_mattype)[m]; }
    inline void set_mattype(const int m, const int mt) { (*_mattype)[m] = mt; }

    inline const double *get_eta(int m) const { return (*_eta)[m]; }
    inline void set_eta( const int i, const double r[NDIMS] );

private:
    // Didn't create a data type for an individual marker to follow the "structure of arrays" concept.
    // Number of markers (may change during remeshing)

    int _nmarkers;
    int _reserved_space;
    int _last_id;

    // Barycentric (local) coordinate within the reference element
    shapefn *_eta;
    // Containing element
    int_vec *_elem;
    // Material type
    int_vec *_mattype;
    // Unique id
    int_vec *_id;

    void random_eta( double* );
    void append_marker( double *eta, int el, int mt );
    void random_markers( const Param&, Variables& );
    void allocate_markerdata( const int );

};

void remap_markers(const Param&, Variables &, 
                   const array_t &, const conn_t &);

#endif
