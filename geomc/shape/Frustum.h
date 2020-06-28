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
 * @brief An N-dimensional frustum (truncated pyramid).
 * 
 * The first N-1 dimensions have cross-sections which are rectangles. The last
 * dimension is considered the "height" of the frustum. 
 * 
 * Frustums need not be symmetrical about the origin. The extents of the cross
 * section are described by an N-1 rectangle lying on the `height = 1` plane.
 *  
 * If the height range spans `height = 0`, the minimum height is effectively 
 * clamped to 0.
 */
template <typename Shape>
class Frustum : public virtual Convex<typename Shape::T, Shape::N + 1> {
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
         * Frustum-point intersection test.
         * 
         * @param p A point.
         * @return `true` if `p` is on or inside this frustum; `false` otherwise.
         */
        bool contains(Vec<T,N> p) const {
            if (not height.contains(p[N-1])) return false;
            p /= p[N-1];
            return base.contains(p.template resized<N-1>());
        }
        
        
        Vec<T,N> convex_support(Vec<T,N> d) const {
            typedef PointType<T,N-1> Pt;
            Rect<T,1> h = clipped_height();
            
            // find the appropriate direction to check in shape-space
            Vec<T,N-1> d_s = d.template resized<N-1>();
            if      (d_s.isZero()) d_s[0] =  1;   // ← pick an arbitrary direction
            else if (h.lo < 0)     d_s    = -d_s; // ← shape is flipped when below origin
            
            // find the relevant point on the base shape, and place it on the h=1 plane
            Vec<T,N> p(base.convex_support(d_s), 1);
            
            // top or bottom face?
            T w = p.dot(d) > 0 ? h.hi : h.lo;
            
            // rescale the point appropriately:
            return w * p;
        }
        
        Rect<T,N> bounds() {
            return base.bounds() * clipped_height();
        }
        
        /// Return the height range of this Frustum after it has been clipped by the origin.
        inline Rect<T,1> clipped_height() const {
            const Rect<T,1>& h = height; // shorthand
            return Rect<T,1>((h.lo < 0 and h.hi > 1) ? 0 : h.lo, h.hi);
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

/// Convenience typedef for oriented, rectangular, N-dimensional Frustums
template <typename T, index_t N>
using ViewFrustum = Oriented< Frustum< Rect<T,N-1> > >;


} // namespace geom

#endif	/* FRUSTUM_H */

