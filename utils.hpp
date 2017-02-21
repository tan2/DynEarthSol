#ifndef DYNEARTHSOL3D_UTILS_HPP
#define DYNEARTHSOL3D_UTILS_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "parameters.hpp"

/////////////////////////////////////////////////////////////////////

class ElemFunc  // base class for functor used in loop_all_elem()
{
public:
    virtual void operator()(int e) = 0;
    virtual ~ElemFunc() {};
};


inline void loop_all_elem(const std::vector<int> &egroups, ElemFunc &functor)
{
#ifdef USE_OMP
    // See mesh.cxx::create_elem_groups() for parallel strategy

    // loop over elements in even element groups
    #pragma omp parallel for default(none) shared(egroups, functor)
    for (std::size_t i=0; i<egroups.size()-1; i+=2) {
        for (int e=egroups[i]; e<egroups[i+1]; ++e)
            functor(e);
    }

    // loop over elements in odd element groups
    #pragma omp parallel for default(none) shared(egroups, functor)
    for (std::size_t i=1; i<egroups.size()-1; i+=2) {
        for (int e=egroups[i]; e<egroups[i+1]; ++e)
            functor(e);
    }

#else

    // loop over all elements sequentially
    for (int e=egroups[0]; e<egroups[egroups.size()-1]; ++e)
        functor(e);

#endif
}

/////////////////////////////////////////////////////////////////////


static void print(std::ostream& os, const double& x)
{
  os << x;
}


static void print(std::ostream& os, const int& x)
{
  os << x;
}


static void print(std::ostream& os, const std::size_t& x)
{
  os << x;
}


template <typename T1, typename T2>
void print(std::ostream& os, const std::pair<T1,T2>& x)
{
    os << x.first << ':' << x.second;
}


template <typename T>
void print(std::ostream& os, const T& A, std::size_t size)
{
  os << "[";
  for (std::size_t i = 0; i != size; ++i) {
    print(os, A[i]);
    if (i+1 != size)
      os << ", ";
  }
  os << "]";
}


template <typename Array>
void print(std::ostream& os, const Array& A)
{
  typename Array::const_iterator i;
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
    os << ", ";
  }
  os << "]";
}


static double trace(const double* s)
{
#ifdef THREED
    return s[0] + s[1] + s[2];
#else
    return s[0] + s[1];
#endif
}


static double second_invariant2(const double* t)
{
#ifdef THREED
    double a = (t[0] + t[1] + t[2]) / 3;
    return ( 0.5 * ((t[0]-a)*(t[0]-a) + (t[1]-a)*(t[1]-a) + (t[2]-a)*(t[2]-a))
             + t[3]*t[3] + t[4]*t[4] + t[5]*t[5] );
#else
    return 0.25*(t[0]-t[1])*(t[0]-t[1]) + t[2]*t[2];
#endif
}


static double second_invariant(const double* t)
{
    /* second invariant of the deviatoric part of tensor t
     * defined as: td = deviatoric(t); sqrt( td(i,j) * td(i,j) / 2)
     */
    return std::sqrt(second_invariant2(t));
}

// C++ verstion of msum() given in
// http://code.activestate.com/recipes/393090-binary-floating-point-summation-accurate-to-full-p/
//
// Description from the original source:
// "Full precision summation using multiple floats for intermediate values"
//   Rounded x+y stored in hi with the round-off stored in lo.  Together
//   hi+lo are exactly equal to x+y.  The inner loop applies hi/lo summation
//   to each partial so that the list of partial sums remains exact.
//   Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:
//   www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps
#if 0
static double accurate_sum(std::vector<double> &a)
{
    double xx, yy, msum = 0.;
    std::vector<double> partials;

    for (auto& xx : a) {
		double x = xx;
        int i = 0;
        for(auto& yy : partials) {
	    	double y = yy;
            if (std::fabs(x) < std::fabs(y)) std::swap(x,y);
            double hi = x + y;
            double lo = y - (hi - x);
            //std::cout << "lo=" << lo <<"\n";
            //std::cout << "hi=" << hi <<"\n";
            if (lo) {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        partials.resize(i+1);
        partials[i] = x;
        //std::cout << "i="<< i << " size="<< partials.size()<< " partials=[";
        //for(auto& el: partials) std::cout << el << ' ';
        //std::cout <<"]"<< std::endl;
    }
    for (auto & n : partials)
		msum += n;

    return msum + 0.0;
}
#endif
static double accurate_sum(double_tbb_vec &a)
{
    double xx, yy, msum = 0.;
    std::vector<double> partials;

    for (auto& xx : a) {
        double x = xx;
        int i = 0;
        for(auto& yy : partials) {
            double y = yy;
            if (std::fabs(x) < std::fabs(y)) std::swap(x,y);
            double hi = x + y;
            double lo = y - (hi - x);
            //std::cout << "lo=" << lo <<"\n";
            //std::cout << "hi=" << hi <<"\n";
            if (lo) {
                partials[i] = lo;
                i += 1;
            }
            x = hi;
        }
        partials.resize(i+1);
        partials[i] = x;
        //std::cout << "i="<< i << " size="<< partials.size()<< " partials=[";
        //for(auto& el: partials) std::cout << el << ' ';
        //std::cout <<"]"<< std::endl;
    }
    for (auto & n : partials)
    msum += n;

    return msum + 0.0;
}
#endif
