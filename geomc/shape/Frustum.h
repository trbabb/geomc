/* 
 * File:   Frustum.h
 * Author: tbabb
 *
 * Created on November 8, 2014, 11:09 PM
 */

#ifndef FRUSTUM_H
#define	FRUSTUM_H

#include <climits>
#include <geomc/shape/Rect.h>
#include <geomc/linalg/AffineTransform.h>

namespace geom {

/**
 * @ingroup shape
 * @brief An N-dimensional frustum (truncated pyramid) with an arbitrary
 * Convex shape as its base, and its (possibly excluded) point at the origin.
 * 
 * The first N-1 dimensions have cross-sections which are `Shape`s. The last
 * dimension is the "height" of the frustum. The base `Shape` is specified 
 * lying on the `h = 1` plane. (Note that below the origin, the base shape
 * will be flipped due to the change of sign).
 *  
 * If the height range spans `height = 0`, then the height extending below the 
 * origin is excluded from the shape. Use `clipped_height()` to obtain the 
 * height range actually spanned by the shape. (A frustum spanning the origin
 * would not be convex).
 *
 * Oriented, Rect-based Frustums are very commonly needed to represent
 * viewing frustums; for these consider using `ViewFrustum`, which is a
 * templated type alias for `Oriented<Frustum<Rect>>`. 
 */
template <typename Shape>
class Frustum : public virtual Convex<typename Shape::elem_t, Shape::N + 1> {
    public:
        /// The coordinate type of this Shape
        typedef typename Shape::elem_t T;
        /// The dimension of this Shape
        static constexpr size_t N = Shape::N + 1;
        
        /// Cross-section of this frustum at the `h = 1` plane.
        Shape base;
        /// Height range spanned by this frustum.
        Rect<T,1> height;
        
        
        /**
         * @brief Construct a pyramid with its tip at the origin
         * and a default-constructed cross-sectional base at h=1.
         */
        Frustum():height(0,1) {}
        
        /**
         * @brief Construct a new Frustum having `base` as a cross section,
         * and spanning heights between `h0` and `h1`.
         */
        Frustum(const Shape& base, T h0, T h1):
                    base(base), 
                    height(std::min(h0,h1), std::max(h0,h1)) {}
        
        /**
         * @brief Construct a new Frustum having `base` as a cross section,
         * and spanning the height range `h`.
         */
        Frustum(const Shape& base, const Rect<T,1>& h):
            base(base),
            height(h) {}
        
        /**
         * Frustum-point intersection test.
         * 
         * @param p A point.
         * @return `true` if `p` is on or inside this frustum; `false` otherwise.
         */
        bool contains(Vec<T,N> p) const {
            if (not clipped_height().contains(p[N-1])) return false;
            if (p[N-1] == 0) return p.isZero();
            p /= p[N-1];
            return base.contains(p.template resized<N-1>());
        }
        
        
        Vec<T,N> convex_support(Vec<T,N> d) const {
            Rect<T,1> h = clipped_height();
            T sign = h.lo < 0 ? -1 : 1; // ← the shape is flipped below the origin
            
            // find the appropriate direction to check in shape-space.
            Vec<T,N-1> d_s = sign * d.template resized<N-1>();
            if (d_s.isZero()) d_s[0] = 1;  // ← pick an arbitrary direction
            
            // find the relevant point on the base shape, and place it on the h=1 plane
            Vec<T,N> p(base.convex_support(d_s), 1);
            
            // top or bottom face?
            // the support point definitely lies on the line passing through 
            // the origin and `p`. which of the two vertices on that line is the 
            // point? the one in the +p direction, or the -p direction?
            T z = (p.dot(d) > 0) ? h.hi : h.lo;
            
            // rescale the point appropriately:
            return z * p;
        }
        
        
        Rect<T,N> bounds() const {
            Rect<T,N-1> b0 = base.bounds();
            Rect<T,1> h    = clipped_height();
            // make two rects for the top and bottom faces;
            // union them; extrude by the height:
            return ((b0 * h.lo) | (b0 * h.hi)) * h; // purdy!!
        }
        
        
        /// Return the height range of this Frustum after it has been clipped by the origin.
        inline Rect<T,1> clipped_height() const {
            const Rect<T,1>& h = height; // shorthand
            return Rect<T,1>((h.lo < 0 and h.hi > 0) ? 0 : h.lo, h.hi);
        }

}; // class Frustum

/**
 * @brief Convenience function to raise the shape `s` into a frustum between 
 * heights `h0` and `h1`, by wrapping `s` in the `Frustum` template.
 *
 * @related Frustum
 */
template <typename Shape>
inline Frustum<Shape> frustum(
        const Shape& s,
        typename Shape::elem_t h0,
        typename Shape::elem_t h1)
{
    return Frustum<Shape>(s, h0, h1);
}

/**
 * @brief Convenience typedef for oriented, rectangular, N-dimensional Frustums
 * @related Frustum
 * @related Oriented
 */
template <typename T, index_t N>
using ViewFrustum = Oriented< Frustum< Rect<T,N-1> > >;


} // namespace geom

#endif	/* FRUSTUM_H */

