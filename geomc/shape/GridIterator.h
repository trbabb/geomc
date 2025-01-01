#pragma once

// todo: rewrite this to use an integer state and convert to/from point_t.

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include <geomc/linalg/Vec.h>
#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Rect.h>

namespace geom {

/***********************************************
 * Grid Iterator                               *
 *                                             *
 * Iterate over all integer points within the  *
 * supplied region.                            *
 ***********************************************/

// TODO: somehow add docs for the inherited operators.
// xxx: this is apparently broken for negative regions.
//    (it's probably a truncating division issue)
    
/**
 * @ingroup shape
 * @brief Iterator over the integer points in an N-dimensional grid.
 * 
 * @tparam T Type of coordinate (integer type recommended).
 * @tparam N Dimension of rectangle to iterate over.
 * @tparam Order Row- or column-major iteration.
 * 
 * This class accepts an N-dimensional rectangular region to iterate over, and returns
 * points on the unit-spaced grid within than region. Row- or column-major order 
 * may be selected with the ArrayOrder template parameter.
 * 
 * GridIterators support all the standard operators expected of a random-access 
 * iterator:
 * 
 *     GridIterator<int,3> i = ... ;
 *     Vec<int,3> pt = *i;
 *     ++i; 
 *     --i;
 *     i += 10;
 *     i1 - i2;
 *     // etc.
 */
template <typename T, index_t N, ArrayOrder Order>
class GridIterator : 
    public boost::iterator_facade
    <
        GridIterator<T,N,Order>,                 // self type
        typename PointType<T,N>::point_t,        // value (pointed to) type
        std::random_access_iterator_tag,         // implemented concepts (all)
        const typename PointType<T,N>::point_t&  // reference to elem type (may be a proxy)
    >
{
    
    static constexpr bool RowMaj = (Order == ARRAYORDER_FIRST_DIM_CONSECUTIVE);
    static constexpr index_t dim_first     = RowMaj ?  0      : (N - 1);
    static constexpr index_t dim_last      = RowMaj ? (N - 1) :  0;
    static constexpr index_t dim_end       = RowMaj ?  N      : -1;
    static constexpr index_t dim_increment = RowMaj ?  1      : -1;
    
public:
    /// Array order (row-major or column-major)
    static constexpr ArrayOrder order = Order;
    /// Type of coordinate to be iterated over.
    typedef typename PointType<T,N>::point_t point_t;
    
    /// Region over which to iterate
    const Rect<T,N> region;
    /// Current point
    point_t pt;
    
    /**
     * Construct an iterator over the region `r` pointing to the first grid point.
     */
    GridIterator(const Rect<T,N> &r):
                region(r),
                pt(region.lo) {}
    /**
     * Construct an interator over the region `r` pointing at the point given by `p`.
     * @param r A region
     * @param p A point inside `r`.
     */
    GridIterator(const Rect<T,N> &r, const point_t &p):
                    region(r),
                    pt(p) {}
    
    /**
     * Construct an iterator over the region bounded by the points
     * `lo` and `hi`, pointing at the first cell.
     * @param lo Lower extreme of the region (inclusive).
     * @param hi Upper extreme of the region (exclusive).
     */
    GridIterator(const point_t &lo, const point_t &hi):
                region(lo, hi - point_t(1)),
                pt(lo) {}
    
    /**
     * @return An iterator pointing to the lower-most coordinate of the grid
     */
    inline GridIterator<T,N,Order> begin() const {
        GridIterator<T,N,Order> other = *this;
        other.pt = region.lo;
        return other;
    }
    
    /**
     * @return An iterator pointing just beyond the upper-most coordinate of the
     * grid.
     */
    inline GridIterator<T,N,Order> end() const {
        GridIterator<T,N,Order> other = *this;
        other.pt = region.lo;
        // one beyond the highest cell. i.e. first cell in
        // (N-1)D block just beyond the end of this region. 
        other.pt[dim_last] = region.hi[dim_last] + 1;
        return other;
    }
    
private:
    
    friend class boost::iterator_core_access;
    
    inline bool equal(const GridIterator<T,N,Order> &other) const {
        return pt == other.pt;
    }
    
    inline void increment() { 
        for (index_t i = dim_first; 
                i != dim_end and ++(pt[i]) > region.hi[i]; 
                i += dim_increment) {
            if (i != dim_last) pt[i] = region.lo[i];
        }
    }
    
    inline void decrement() {
        for (index_t i = dim_first;
                i != dim_end and --(pt[i]) < region.lo[i];
                i += dim_increment) {
            if (i != dim_last) pt[i] = region.hi[i] - 1;
        }
    }
    
    inline void advance(index_t n) {
        point_t dest = point_t(0);
        point_t dim  = region.dimensions();
        index_t vol  = region.volume();
        
        for (index_t i = dim_last; i != dim_first - dim_increment; i -= dim_increment) {
            vol     = vol / dim[i];
            dest[i] = n / vol;
            n       = n % vol;
        }
        
        pt = dest;
    }
    
    inline index_t distance_to(const GridIterator<T,N,Order> &other) const {
        point_t dim  = region.dimensions(); 
        point_t jump = other.pt - pt;
        index_t dist = 0;
        index_t vol  = 1;
        for (index_t i = dim_first; i != dim_end; i += dim_increment) {
            dist += vol * jump[i];
            vol  *= std::max(dim[i], (T)0);
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
class GridIterator<T,1,Order> : 
    public boost::iterator_facade
    <
        GridIterator<T,1,Order>,                 // self type
        typename PointType<T,1>::point_t,        // value (pointed to) type
        std::random_access_iterator_tag,         // implemented concepts (all)
        const typename PointType<T,1>::point_t&  // reference to elem type (may be a proxy)
    >
{

public:

    static const ArrayOrder order;
    typedef typename PointType<T,1>::point_t point_t;
    
    const Rect<T,1> region;
    point_t pt;
    
    GridIterator(const Rect<T,1> &r):
            region(r),
            pt(region.lo) {}
    GridIterator(const point_t &lo, const point_t &hi):
            region(lo, hi + 1),
            pt(lo) {}
    
    inline GridIterator<T,1,Order> begin() const {
        GridIterator<T,1> other = *this;
        other.pt = region.lo;
        return other;
    }
    
    inline GridIterator<T,1,Order> end() const {
        GridIterator<T,1,Order> other = *this;
        other.pt = region.hi + 1;
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

} /* namespace geom */
