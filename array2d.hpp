#if !defined(DYNEARTHSOL3D_ARRAY2D_h)
#define DYNEARTHSOL3D_ARRAY2D_h

#include <algorithm>


template <typename T, int N>
class Array2D {

    T* a_;
    int n_;
    // Can this object safely delete a_?
    // When a_ is owned by another object, or when a_ is mallac'd (not
    // returned by the new operator), we cannot delete a_ in this object.
    bool can_delete;

public:
    //
    // constructors & destructor
    //
    Array2D() {a_ = NULL; n_ = 0; can_delete = 1;}
    Array2D(T* a, int n) {a_ = a; n_ = n; can_delete = 0;}

    Array2D(int size, const T& val) {
        a_ = new T[N*size];
        n_ = size;
        std::fill_n(a_, N*n_, val);
        can_delete = 1;
    }

    explicit
    Array2D(int size) {
        a_ = new T[N*size];
        n_ = size;
        can_delete = 1;
    }

    ~Array2D() {if (can_delete) delete [] a_;}

    //
    // methods
    //
    T* data() {return a_;}
    const T* data() const {return a_;}
    std::size_t size() const {return n_;}
    int num_elements() const {return N*n_;}

    // NOTE: reset DOES NOT delete the array
    void reset() {a_ = NULL; n_ = 0; can_delete = 0;}

    // steal the pointer from other
    void steal_ref(Array2D& other) {
        if (can_delete) delete [] a_;
        a_ = other.data();
        n_ = other.size();
        can_delete = other.can_delete;
        other.reset();
    }

    //
    // index accessing
    //
    T* operator[](std::size_t i) {return a_ + N*i;}
    const T* operator[](std::size_t i) const {return a_ + N*i;}

    //
    // iterators
    //
    typedef T* iterator;
    typedef const T* const_iterator;
    iterator begin() {return a_;}
    const_iterator begin() const {return a_;}
    iterator end() {return a_ + N*n_;}
    const_iterator end() const {return a_ + N*n_;}

    typedef T element;

private:
    // disable copy and assignment operators
    Array2D(const Array2D&);
    Array2D<T,N>& operator=(const Array2D<T,N>& rhs);
};

#endif
