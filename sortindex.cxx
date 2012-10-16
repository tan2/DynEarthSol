#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>


class RangeGenerator {
    int current;
public:
    RangeGenerator(int init=0) {current=init;}
    int operator()() {return current++;}
};


template< typename T >
class idx_lt {
    const std::vector<T>& _x;
public:
    idx_lt( const std::vector<T>& x ) : _x(x) {}
    bool operator()( std::size_t left, std::size_t right ) const { return _x[left] < _x[right]; }
};


template< typename T >
void sortindex(const std::vector<T>& x, std::vector<std::size_t>& idx)
{
    // fill idx = [0, 1, 2, ...]
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), idx_lt<T>(x));
}


void test_range()
{
    RangeGenerator range;
    for(int i=0; i<10; ++i) {
        std::cout << range() << '\n';
    }

    RangeGenerator range1(100);
    for(int i=0; i<10; ++i) {
        std::cout << range1() << '\n';
    }
}


void test_sort()
{
    RangeGenerator range1(3), range2;
    std::vector<std::size_t> list(10), idx(10);
    std::generate(list.begin(), list.begin()+5, range1);
    std::generate(list.begin()+5, list.end(), range2);

    sortindex(list, idx);

    std::cout << "result ...\n";
    for(int i=0; i<10; ++i) {
        std::cout << list[i] << ", " << idx[i] << '\n';
    }

    std::cout << "sorted ...\n";
    for(int i=0; i<10; ++i) {
        std::cout << i << ": " << list[idx[i]] << '\n';
    }
}
