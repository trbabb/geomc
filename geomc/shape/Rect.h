/*
 * Rect.h
 * 
 *  Created on: Oct 7, 2010
 *      Author: Tim Babb
 */


#ifndef Rect_H_
#define Rect_H_

#include <algorithm>
#include <cmath>
#include <limits>

#if __cplusplus < 201103L
#include <tr1/functional>
#else
#include <functional>
#endif

#include <boost/utility/enable_if.hpp>

#include <geomc/Hash.h>
#include <geomc/shape/Bounded.h>
#include <geomc/linalg/Ray.h>
#include <geomc/linalg/Vec.h>
#include <geomc/shape/shapedetail/Hit.h>

// all boundaries are inclusive. this is inconsistent with half-open interval convention, but:
// - union of rect with a point should thereafter contain the point
// - with a half-open interval, we'd need to "increment" the upper bound
//   to contain p
// - which is slow and ugly
// - the upper and lower bounds now mean different things, and
//   that means you need to increment and decrement
//   Sometimes And Not Other Times. 
// - compare to: "both endpoints are symmetrical and the same logic applies everywhere"
// - the biggest downside: ownership of points on the faces of abutting Rects
//   is ambiguous (for Rects with real-valued coordinates).

namespace geom {

/****************************
 * Rect Class               *
 ****************************/

/** @ingroup shape
 *  @brief An N-dimensional axis-aligned interval.
 *
 * This class works naturally with both 1D and multi-dimensional ranges.
 * 1D ranges have point type `T`, while multi-dimensional ranges have point
 * type `Vec<T,N>`. A typedef of this may be accessed via:
 * 
 *     Rect<T,N>::point_t
 *
 * Example: 
 * 
 * * `Rect<double,4>::point_t` is `Vec<double,4>`; `lo` and `hi` are `Vec4d`.
 * * `Rect<double,1>::point_t` is `double`; `lo` and `hi` are `double`.
 * 
 * Rect boundaries are inclusive. The points `lo` and `hi` are considered
 * to be inside the Rect. Therefore, "degenerate" Rects (like those that contain only
 * a single point, edge, or face), are considered non-empty.
 *
 * Note that this convention differs slightly from conventional "interval" logic
 * on integers, wherein the upper boundary is excluded. This convention was chosen
 * to make all the faces of a Rect symmetrical, and avoid special cases or conditional
 * increments. To iterate over the range within an integer Rect, consider using a GridIterator
 * object, which abstracts away the boundary logic for you.
 *
 * If any coordinate of `lo` is greater than the same coordinate of `hi`, then
 * the `Rect` is empty.
 */
template <typename T, index_t N> 
class Rect : virtual public Convex<T,N> {
protected:
    typedef PointType<T,N> ptype;

public:
    /// Type of object confined by this Rect
    typedef typename PointType<T,N>::point_t point_t;
    
    /// Lower extremes
    point_t lo;
    /// Upper extremes
    point_t hi;
    
    static const point_t endpoint_measure;

public:

    /**
     * @brief Construct a Rect with extremes `lo` and `hi`. 
     * If for any axis `lo > hi`, the Rect is empty.
     * 
     * @param lo Lower extreme
     * @param hi Upper extreme
     */
    Rect(point_t lo, point_t hi):
        lo(lo),
        hi(hi) {}
    
    
    /**
     * @brief Construct a Rect containing only the point `p`.
     */
    Rect(point_t p):
        lo(p),
        hi(p) {}
    
    
    /**
     * @brief Construct an empty interval.
     *
     * Sets lower and upper bounds to the maximum and minimum values of `T`, respectively.
     *
     * A union between this Rect and any finite shape is an identity operation.
     */
    Rect():
        lo(std::numeric_limits<T>::max()),
#if __cplusplus >= 201103L
        hi(std::numeric_limits<T>::lowest()) {
#else
        // c++03, you make me sad.
        hi(  std::numeric_limits<T>::is_integer ? 
               std::numeric_limits<T>::min() : 
              -std::numeric_limits<T>::max()) {
#endif
        //do nothing else
    }


    /****************************
     * Convenience functions    *
     ****************************/

    /**
     * @brief Construct a Rect from a center point and extent.
     * @param c Center of the new Rect
     * @param dims Lengths of each axis; may be negative.
     * @return A new Rect with its center at `c`.
     */
    inline static Rect<T,N> from_center(
            typename Rect<T,N>::point_t c,
            typename Rect<T,N>::point_t dims)
    {
        return Rect<T,N>::spanning_corners(
            c -  dims / 2, 
            c + (dims - endpoint_measure) / 2);
    }
    
    /**
     * @brief Construct a Rect from a corner point and an extent.
     * @param corner Arbitrary corner point.
     * @param dims Lengths of each axis relative to given corner. Lengths may be negative.
     * @return A new Rect with one corner at `corner`.
     */
    inline static Rect<T,N> from_corner(
            typename Rect<T,N>::point_t corner,
            typename Rect<T,N>::point_t dims)
    {
        return Rect<T,N>::spanning_corners(
            corner, 
            corner + dims - endpoint_measure);
    }
    
    /**
     * @brief Construct a Rect containing the two corners `c1` and `c2`. 
     * @param c1 A corner of the Rect 
     * @param c2 The corner opposite `c1`
     * @return A new Rect with corners at `c1` and `c2`.
     */
    inline static Rect<T,N> spanning_corners(
            point_t c1,
            point_t c2)
    {
        typename Rect<T,N>::point_t lo = std::min(c1, c2);
        typename Rect<T,N>::point_t hi = std::max(c1, c2);
        return Rect<T,N>(lo, hi);
    }
    
    /**
     * @brief Construct a Rect containing all the points in the sequence
     * between `begin` and `end`.
     * @tparam PointIterator An iterator to a `point_t`.
     * @param begin Iterator to the first point.
     * @param end Iterator to the off-end point.
     * @return A new Rect containing all the points in the given sequence.
     */
    template <typename PointIterator>
    inline static Rect<T,N> from_points(PointIterator begin, PointIterator end) {
        Rect<T,N> r;
        for (PointIterator i = begin; i != end; ++i) {
            r |= *i;
        }
        return r;
    }

    /**
     * @brief Test whether a point is in the N-dimensional range `[lo, hi]`.
     *
     * If `lo` > `hi` along any axis, then the range is empty and the 
     * function returns `false`.
     * 
     * @param lo Lower extreme
     * @param hi Upper extreme
     * @param pt Test point
     * @return `true` if `pt` is inside the exremes
     */
    static bool contains(typename Rect<T,N>::point_t lo,
                         typename Rect<T,N>::point_t hi,
                         typename Rect<T,N>::point_t pt) {
        for (index_t axis = 0; axis < N; axis++) {
            T v = ptype::iterator(pt)[axis];
            if (v < ptype::iterator(lo)[axis] or v > ptype::iterator(hi)[axis]) {
                return false;
            }
        }
        return true;
    }

    /*****************************
     * Operators                 *
     *****************************/

    /**
     * @brief Interval union.
     * 
     * @return A Rect fully containing `this` and box `b`.
     */
    inline Rect<T,N> operator|(const Rect<T,N> &b) const {
        return range_union(b);
    }
    
    /**
     * @brief Interval union.
     * 
     * Extend `this` to fully contain `b`.
     * @return A reference to `this`, for convenience.
     */
    Rect<T,N>& operator|=(const Rect<T,N> &b) {
        //box union
        hi = std::max(b.hi, hi);
        lo = std::min(b.lo, lo);
        return *this;
    }
    
    /**
     * @brief Point union.
     *
     * @return A Rect fully containing `this` and the point `p`.
     */
    Rect<T,N> operator|(const point_t& p) {
        // point union
        return Rect<T,N>(
            std::min(lo, p),
            std::max(hi, p)
        );
    }
    
    /**
     * @brief Point union.
     *
     * Extend `this` to fully contain `p`. 
     * @return a reference to `this`, for convenience.
     */
    Rect<T,N>& operator|=(const point_t& p) {
        // point union
        lo = std::min(lo, p);
        hi = std::max(hi, p);
        return *this;
    }
    
    /**
     * @brief Interval intersection.
     * 
     * @return A Rect representing the area overlapped by both `this` and `b`.
     */
    inline Rect<T,N> operator&(const Rect<T,N> &b) const {
        return range_intersection();
    }
    
    /**
     * @brief Interval intersection.
     * 
     * Confine this Rect to the area overlapping `b`.
     * @return A reference to `this` for convenience.
     */
    Rect<T,N>& operator&=(const Rect<T,N> &b) {
        lo = std::max(lo, b.lo);
        hi = std::min(hi, b.hi);
    }

    /**
     * @brief Translation.
     * 
     * @param dx Amount by which to translate this region.
     *   Add `dx` to the coordinates of all the bounds.
     * @return A translated Rect.
     */
    inline Rect<T,N> operator+(point_t dx) const {
        return Rect<T,N>(lo + dx, hi + dx);
    }
    
    /**
     * @brief Translation.
     * 
     * @param dx Amount by which to translate this region.
     *   Subtract `dx` from the coordinates of all the bounds.
     * @return A translated Rect.
     */
    inline Rect<T,N> operator-(point_t dx) const {
        return Rect<T,N>(lo - dx, hi - dx);
    }
    

    /**
     * @brief Translation.
     * 
     * @param dx Amount by which to translate this region.
     *   Add `dx` to the coordinates of all the bounds.
     * @return  A reference to `this`, for convenience.
     */
    inline Rect<T,N>& operator+=(point_t dx) {
        hi += dx;
        lo += dx;
        return *this;
    }
    
    
    /**
     * @brief Translation.
     * 
     * @param dx Amount by which to translate this region. 
     *   Subtract `dx` from the coordinates of all the bounds.
     * @return  A reference to `this`, for convenience.
     */
    inline Rect<T,N>& operator-=(point_t dx) {
        hi -= dx;
        lo -= dx;
        return *this;
    }
    
    
    /**
     * @brief Interval exclusion.
     *
     * @param other Range to be excluded from the range of this `Rect`.
     * @return A new `Rect` whose range does not overlap `other`.
     */
    inline Rect<T,N> operator-(const Rect<T,N>& other) {
        return Rect<T,N>(
            std::max(lo,other.hi),
            std::min(hi,other.lo));
    }
    
    
    /**
     * @brief Interval exclusion.
     *
     * @param other Range to remove from this `Rect`.
     * @return A reference to `this`, for convenience.
     */
    inline Rect<T,N>& operator-=(const Rect<T,N>& other) {
        lo = std::max(lo, other.hi);
        hi = std::min(hi, other.lo);
        return *this;
    }
    

    /**
     * @brief Equality test.
     * 
     * @return `true` if and only if all the corresponding extremes of `b` are the same.
     */
    inline bool operator==(Rect<T,N> b) const {
        return (hi == b.hi) and (lo == b.lo);
    }
    
    /**
     * @brief Inequality test.
     * 
     * @return `true` if and only if any extreme of `b` is different from
     * the corresponding extreme in `this`.
     */
    inline bool operator!=(Rect<T,N> b) const {
        return (hi != b.hi) or (lo != b.lo);
    }

    /**
     * @brief Scale transformation.
     * 
     * @param a Scale factor.
     * @return A new Rect, scaled about the origin by factor `a`.
     */
    inline Rect<T,N> operator*(point_t a) const {
        return Rect<T,N>(lo * a, hi * a);
    }

    /**
     * @brief Scale transformation. 
     * 
     * Scale this Rect about the origin by factor `a`. 
     * 
     * @param a Scale factor
     * @return A reference to `this`, for convenience.
     */
    inline Rect<T,N>& operator*=(point_t a) {
        lo *= a;
        hi *= a;
        return *this;
    }

    /**
     * @brief Scale transformation.
     * 
     * @param a Scale factor.
     * @return A new Rect, scaled about the origin by multiple `1 / a`.
     */
    inline Rect<T,N> operator/(point_t a) const {
        return Rect<T,N>(lo / a, hi / a);
    }

    /**
     * @brief Scale transformation. 
     * 
     * Scale this Rect about the origin by factor `1 / a`. 
     * @param a Scale factor
     * @return A reference to `this`, for convenience.
     */
    inline Rect<T,N>& operator/=(point_t a) {
        lo /= a;
        hi /= a;
        return *this;
    }
    
    /**
     * @brief Interval cartesian product. 
     *
     * Extrude this Rect into a higher dimension, by treating the extents of 
     * `r` as the extents along the new dimensions. In other words, concatenate 
     * the coordinates of `this` and `r` into a new Rect.
     *
     * For example, the product of a 2D rectangle with a 1D range is a 3D box.
     * If multiplied in that order, the product will have the extents of the rectangle 
     * along its first two axes and the extents of the 1D range along its third axis.
     *
     * @tparam M Dimensionality of `r`.
     */
    template <index_t M>
    inline Rect<T,M + N> operator*(const Rect<T,M>& r) const {
        Rect<T,M+N> o;
        
        std::copy(ptype::iterator(lo),   ptype::iterator(lo)   + N, o.lo.begin());
        std::copy(ptype::iterator(hi),   ptype::iterator(hi)   + N, o.hi.begin());
        std::copy(ptype::iterator(r.lo), ptype::iterator(r.lo) + M, o.lo.begin() + N);
        std::copy(ptype::iterator(r.lo), ptype::iterator(r.hi) + M, o.hi.begin() + N);
        
        return o;
    }

    /**
     * @brief Element-wise typecast.
     * 
     * @return A new Rect, with elements all of type `U`.
     */
    template <typename U, index_t M> operator Rect<U,M>() const {
        return Rect<U,M>((typename Rect<U,M>::point_t) lo,
                         (typename Rect<U,M>::point_t) hi);
    }


    /*****************************
     * Public Methods            *
     *****************************/

    /**
     * @return `true` if and only if `pt` is inside this rectangle.
     * Points on the surface of the Rect are considered to be contained by it.
     */
    bool contains(point_t pt) const {
        for (index_t axis = 0; axis < N; axis++) {
            T v = ptype::iterator(pt)[axis];
            if (v < ptype::iterator(lo)[axis] or 
                v > ptype::iterator(hi)[axis])
            {
                return false;
            }
        }
        return true;
    }
 
    /**
     * @return `true` if and only if there is a point overlapped by both `Rect`s.
     */
    bool intersects(const Rect<T,N> &box) const {
        for (index_t axis = 0; axis < N; axis++) {
            // disjoint on this axis?
            if (ptype::iterator(hi)[axis]     <= ptype::iterator(box.lo)[axis] or 
                ptype::iterator(box.hi)[axis] <= ptype::iterator(lo)[axis]) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * @return The center point of this region.
     */
    point_t center() const {
        return (hi + lo + endpoint_measure) / 2;
    }

    /**
     * @return The size of this region along each axis.
     *
     * Note that for integer type Rects, since both the high and low boundaries are included,
     * the length along each axis is `hi - lo + 1`.
     */
    inline point_t dimensions() const {
        return hi - lo + endpoint_measure;
    }

    /**
     * Change the size of this region, adjusting about its center.
     * @param dim New lengths along each axis.
     */
    void set_dimensions(point_t dim) {
        // rule: preserve the center point.
        //       for odd-length boxes, this is the center element.
        //       for even-length boxes, it is the middle zero-indexed "fencepost".
        // [0|1|2|3]
        //   [1|2|3]
        
        // [0|1|2|3]
        //   [1|2]
        
        // [0|1|2|3|4]
        //   [1|2|3]
        
        // [0|1|2|3|4]
        // [0|1|2|3]
        
        // xxx fixme for int types
        dim = std::abs(dim);
        point_t diff = (dim - dimensions()) / 2;
        lo = lo - diff;
        hi = hi + diff;
    }

    /**
     * Re-configure this region to exactly contain the two given points.
     */
    void set_corners(point_t corner1, point_t corner2) {
        hi = std::max(corner1, corner2);
        lo = std::min(corner1, corner2);
    }
    
    /**
     * Preserving its size, translate the center of this region to the point
     * given by `center`.
     * @param c New center point.
     */
    void center_on(point_t c) {
        point_t tx = c - center();
        hi += tx;
        lo += tx;
    }

    /**
     * Translation
     * @param tx Amount by which to translate this region.
     */
    void translate(point_t tx) {
        hi += tx;
        lo += tx;
    }
    
    /**
     * Compute the N-dimensional volume of this region as the product of its
     * extents. Result measures a length if `N` is 1; an area if `N` is 2; a 
     * volume if 3; etc. 
     * 
     * @return The volume of this region.
     */
    T volume() const {
        T vol = 1;
        point_t dim = dimensions();
        for (index_t axis = 0; axis < N; axis++) {
            vol *= std::max(ptype::iterator(dim)[axis], (T)0);
        }
        return vol;
    }
    
    /**
     * @return `true` if and only if this region contains no points; which
     * will be the case if `lo[i] > hi[i]` for any axis `i`.
     */
    bool is_empty() const {
        for (index_t axis = 0; axis < N; axis++) {
            T hi_i = ptype::iterator(hi)[axis];
            T lo_i = ptype::iterator(lo)[axis];
            if (hi_i < lo_i) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * @return A new Rect representing the area overlapped by both `this` and `b`.
     */
    Rect<T,N> range_intersection(const Rect<T,N> &b) const {
        return Rect<T,N>(std::max(lo, b.lo),
                         std::min(hi, b.hi));
    }

    /**
     * @return A new Rect completely containing both `this` and `b`.
     */
    Rect<T,N> range_union(const Rect<T,N> &b) const {
        return Rect<T,N>(std::min(lo, b.lo),
                         std::max(hi, b.hi));
    }
    
    /**
     * Clamp the coordinates of `p` to lie within this Rect.
     *
     * Result can be considered the point nearest to `p` contained in this `Rect`.
     */
    inline point_t clamp(point_t p) const {
        return std::min(hi, std::max(lo, p));
    }
    
    
    /**
     * @return The square of the distance to the nearest point contained by this `Rect`; 
     * zero if `p` is inside.
     */
    inline T dist2(point_t p) const {
        return ptype::mag2(p - clamp(p));
    }
    
    /**
     * Remap `s` on the `[0,1]` interval to the extents of this Rect. 
     * In other words, interpolate the corners of this Rect using s as an interpolation 
     * parameter.
     *
     * Inverse operation of `unmap()`.
     * 
     * Values of `s` between 0 and 1 correspond to points inside this Rect.
     */
    point_t remap(point_t s) const {
        // todo: template for integer type? beware endpoint measure.
       return (point_t(1) - s) * lo + s * hi;
    }
    
    /**
     * Find `p`'s fractional position within this Rect.
     *
     * Inverse operation of `remap()`.
     */
    point_t unmap(point_t p) const {
        return (p - lo) / (hi - lo);
    }
    
    point_t convex_support(point_t d) const {
        point_t o;
        for (index_t i = 0; i < N; i++) {
            T a = ptype::iterator(d)[i];
            ptype::iterator(o)[i] = ptype::iterator(a < 0 ? lo : hi)[i];
        }
        return o;
    }
    
    /**
     * Box-ray intersection test.
     * 
     * @param r A ray
     * @param sides Whether to hit-test the front-facing or back-facing surfaces.
     * @return A ray Hit describing whether and where the ray intersects this Rect,
     * as well as the normal, side hit, and ray parameter.
     */
    Hit<T,N> trace(const Ray<T,N> &r, HitSide sides) const {
        // use non-denormal extremes, for speed (as opposed to inf).
        // (also, infinity may not be defined for T).
        T near_root = std::numeric_limits<T>::lowest();
        T far_root  = std::numeric_limits<T>::max();
        index_t near_axis = 0;
        index_t far_axis  = 0;
        // for guaranteeing that the hit point is exactly on the rect surface,
        // regardless of precision and rounding error in the ray arithmetic:
        T near_coord  = 0;
        T far_coord   = 0;
        Hit<T,N> hit  = Hit<T,N>(r, sides); // defaults to miss
        
        // test ray intersection with an infinite slab along each axis
        // box intersections are the farthest near hit and nearest far hit
        for (index_t axis = 0; axis < N; axis++) {
            if (r.direction[axis] == 0) {
                // ray direction tangent to test planes, no intersection along this axis
                if (r.origin[axis] < lo[axis] || r.origin[axis] > hi[axis]) {
                    // origin outside of test slab; miss
                    return hit;
                }
            } else {
                // coordinate of hit, along tested axis
                T c1 = hi[axis];
                T c2 = lo[axis];
                // ray multiple of hit
                T s1 = (c1 - r.origin[axis]) / r.direction[axis];
                T s2 = (c2 - r.origin[axis]) / r.direction[axis];
                
                // order s1, s2 to (near, far) (aka front, back)
                if (s2 < s1) {
                    std::swap(s1,s2);
                    std::swap(c1,c2);
                }
                if (s1 > near_root) {
                    near_root  = s1;
                    near_axis  = axis; 
                    near_coord = c1;
                }
                if (s2 < far_root) {
                    far_root  = s2;
                    far_axis  = axis;
                    far_coord = c2;
                }
            }
        }
        
        if (near_root > far_root) {
            // miss
            return hit;
        } else if (near_root > 0 and (sides & HIT_FRONT)) {
            // hit front
            hit.p = r * near_root;
            hit.p[near_axis] = near_coord;
            hit.s = near_root;
            hit.n[near_axis] = -copysign(1, r.direction[near_axis]);
            hit.side = HIT_FRONT;
            hit.hit = true;
        } else if (far_root > 0 and (sides & HIT_BACK)) {
            // hit back
            hit.p = r * far_root;
            hit.p[far_axis] = far_coord;
            hit.s = far_root;
            hit.n[far_axis] = copysign(1, r.direction[far_axis]);
            hit.side = HIT_BACK;
            hit.hit = true;
        }
        // if encountered cases above; hit.
        // else miss (box entirely behind ray origin)
        return hit;
    }

    /*****************************
     * Inherited Methods         *
     *****************************/

    Rect<T,N> bounds() const {
        return *this;
    }

}; // end Rect class



// static members
template <typename T, index_t N>
const typename Rect<T,N>::point_t Rect<T,N>::endpoint_measure = 
    std::numeric_limits<T>::is_integer ? 1 : (T)0;


} // end namespace geom

// allows RectBound<T,N>s to work with std::tr1 hash containers
namespace std {

#if __cplusplus < 201103L
namespace tr1 {
#endif

   template <typename T, index_t N>
   struct hash< geom::Rect<T,N> > : public unary_function<geom::Rect<T,N>, size_t> {
       inline size_t operator() (const geom::Rect<T,N>& box) const {
           return geom::general_hash(&box, sizeof(geom::Rect<T,N>));
       }
   };

#if __cplusplus < 201103L
} // namespace tr1
#endif

} // namespace std

#endif /* Rect_H_ */
