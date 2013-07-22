/*
 * GridIterator.h
 *
 *  Created on: Mar 14, 2013
 *      Author: tbabb
 */

#ifndef GRIDITERATOR_H_
#define GRIDITERATOR_H_

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/static_assert.hpp>
#include "linalg/Vec.h"
#include "shape/ShapeTypes.h"
#include "shape/Rect.h"
#include "shape/shapedetail/GridDetail.h"

#define ENABLEINT_TYPE(T) boost::enable_if<boost::is_integral<T>, int>::type

namespace geom {

/***********************************************
 * Grid Iterator                               *
 *                                             *
 * Iterate over all integer points within the  *
 * supplied region.                            *
 ***********************************************/

// TODO: should this really be templated over array order? maybe that should be dynamic?
// TODO: This shit is broken and walking off the edge of things (reverse dimension order only).

template <typename T, index_t N, ArrayOrder Order>
class GridIterator : public boost::iterator_facade<GridIterator<T,N,Order>,      // self type
                                            typename PointType<T,N>::point_t,          // value (pointed to) type
                                            std::random_access_iterator_tag,           // implemented concepts (all)
                                            const typename PointType<T,N>::point_t&> { // reference to elem type (may be a proxy)
    
    const static index_t dim_first     = detail::_ImplArrayOrder<Order,N>::dim_first;
    const static index_t dim_last      = detail::_ImplArrayOrder<Order,N>::dim_last;
    const static index_t dim_end       = detail::_ImplArrayOrder<Order,N>::dim_end;
    const static index_t dim_increment = detail::_ImplArrayOrder<Order,N>::dim_increment;
    
public:
    const static ArrayOrder order; //for cxx11: use constexpr
    typedef typename PointType<T,N>::point_t point_t;
    
    const Rect<T,N> region;
    point_t pt;
    
    GridIterator(const Rect<T,N> &r):
                region(r),
                pt(region.min()) {}
    GridIterator(const Rect<T,N> &r, const point_t &p):
                    region(r),
                    pt(p) {}
    GridIterator(const point_t &lo, const point_t &hi):
                region(lo,hi),
                pt(lo) {}
    
    inline GridIterator<T,N,Order> begin() const {
        GridIterator<T,N,Order> other = *this;
        other.pt = region.min();
        return other;
    }
    
    inline GridIterator<T,N,Order> end() const {
        GridIterator<T,N,Order> other = *this;
        other.pt = region.min();
        // one beyond the highest cell. i.e. first cell in
        // (N-1)D block just beyond the end of this region. 
        other.pt[dim_last] = region.max()[dim_last];
        return other;
    }
    
private:
    
    friend class boost::iterator_core_access;
    
    inline bool equal(const GridIterator<T,N,Order> &other) const {
        return pt == other.pt;
    }
    
    inline void increment() { 
        for (index_t i = dim_first; i != dim_end; i += dim_increment){
            if (pt[i] >= region.max()[i] - 1) {
                pt[i] = region.min()[i];
            } else {
                pt[i]++;
                return;
            }
        }
        // just beyond end.
        pt[dim_last] = region.max()[dim_last];
    }
    
    inline void decrement() {
        for (index_t i = dim_first; i != dim_end; i += dim_increment){
            if (pt[i] <= region.min()[i]) {
                pt[i] = region.max()[i] - 1;
            } else {
                pt[i]--;
                return;
            }
        }
        // if we don't do this, then we'll wrap around again to max();
        // we can never leave the region. i == i.end() needs to be
        // possible, so we allow the pointer to increment/decrement 
        // beyond the edge.
        pt[dim_last] = region.min()[dim_last] - 1;
    }
    
    inline void advance(index_t n) {
        point_t dest = point_t(0);
        point_t dim  = region.getDimensions();
        index_t vol  = region.volume();
        
        for (index_t i = dim_last; i != dim_first - dim_increment; i -= dim_increment) {
            vol     = vol / dim[i];
            dest[i] = n / vol;
            n       = n % vol;
        }
        
        pt = dest;
    }
    
    inline index_t distance_to(const GridIterator<T,N,Order> &other) const {
        // TODO: test this. handling of volume's sign is suspicious.
        // TODO: handling of reverse case is suspicious too.
        point_t dim  = region.getDimensions(); 
        point_t jump = other.pt - pt;
        index_t dist = 0;
        index_t vol  = 1;
        for (index_t i = dim_first; i != dim_end; i += dim_increment) {
            dist += vol * jump[i];
            vol  *= dim[i];
        }
        return dist;
    }
    
    inline const point_t& dereference() const {
        return pt;
    }
    
};


/***********************************************
 * N = 1 specialization                        *
 ***********************************************/

// a pretty stupid class, but ensures that other templated classes
// that use a GridIterator still work as expected.

// probably not interoperable with other orderings, though the
// classes are identical. Could be fixed, perhaps, with a templated 
// conversion operator or constructor.

template <typename T, ArrayOrder Order>
class GridIterator<T,1,Order> : public boost::iterator_facade<GridIterator<T,1,Order>,      // self type
                                                 typename PointType<T,1>::point_t,          // value (pointed to) type
                                                 std::random_access_iterator_tag,           // implemented concepts (all)
                                                 const typename PointType<T,1>::point_t&> {    // reference to elem type (may be a proxy)
    
public:

    static const ArrayOrder order;
    typedef typename PointType<T,1>::point_t point_t;
    
    const Rect<T,1> region;
    point_t pt;
    
    GridIterator(const Rect<T,1> &r):
            region(r),
            pt(region.min()) {}
    GridIterator(const point_t &lo, const point_t &hi):
            region(lo,hi),
            pt(lo) {}
    
    inline GridIterator<T,1,Order> begin() const {
        GridIterator<T,1> other = *this;
        other.pt = region.min();
        return other;
    }
    
    inline GridIterator<T,1,Order> end() const {
        GridIterator<T,1,Order> other = *this;
        other.pt = region.max();
        return other;
    }
    
private:
    
    friend class boost::iterator_core_access;
    
    inline bool equal(const GridIterator<T,1,Order> &other) const {
        return pt == other.pt;
    }
    
    inline void increment() { 
        pt++;
    }
    
    inline void decrement() {
        pt--;
    }
    
    inline void advance(index_t n) {
        pt += n;
    }
    
    inline index_t distance_to(const GridIterator<T,1,Order> &other) const {
        return other.pt - pt;
    }
    
    inline const point_t& dereference() const {
        return pt;
    }
    
};

// static variable for all specializations

template <typename T, index_t N, ArrayOrder Order>
const ArrayOrder GridIterator<T,N,Order>::order = Order;

}; /* namespace geom */

#undef ENABLEINT_TYPE

#endif /* GRIDITERATOR_H_ */
