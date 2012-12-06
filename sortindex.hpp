#ifndef DYNEARTHSOL3D_SORTINDEX_HPP
#define DYNEARTHSOL3D_SORTINDEX_HPP

#include <numeric>

template< typename T >
class idx_lt {
    const std::vector<T>& _x;
public:
    idx_lt( const std::vector<T>& x ) : _x(x) {}
    bool operator()( std::size_t left, std::size_t right ) const { return _x[left] < _x[right]; }
};


template< typename T >
class idx_gt {
    const std::vector<T>& _x;
public:
    idx_gt( const std::vector<T>& x ) : _x(x) {}
    bool operator()( std::size_t left, std::size_t right ) const { return _x[left] > _x[right]; }
};


/* sort x in ascending order, x is not modified, the sorted order is stored in idx */
template< typename T >
void sortindex(const std::vector<T>& x, std::vector<std::size_t>& idx)
{
    // fill idx = [0, 1, 2, ...]
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), idx_lt<T>(x));
}


/* sort x in descending order, x is not modified, the sorted order is stored in idx */
template< typename T >
void sortindex_reversed(const std::vector<T>& x, std::vector<std::size_t>& idx)
{
    // fill idx = [0, 1, 2, ...]
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), idx_gt<T>(x));
}


#endif
