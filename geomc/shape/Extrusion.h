/* 
 * File:   Frustum.h
 * Author: tbabb
 *
 * Created on November 8, 2014, 11:09 PM
 */

#ifndef EXTRUSION_H
#define	EXTRUSION_H

#include <climits>
#include <geomc/shape/Rect.h>
#include <geomc/linalg/AffineTransform.h>


namespace geom {

/**
 * @ingroup shape
 * @brief An N-dimensional extrusion of an arbitrary N-1 dimensional Convex shape.
 * 
 * The first N-1 dimensions have cross-sections which are `Shape`s
 * (e.g., `Rect`, `Sphere`, `Frustum`...). The shape is lofted "vertically" into 
 * the last dimension.
 *
 * Example:
 * 
 *     // initialize a cylinder of length and radius 1:
 *     auto cylinder = Extrusion<Sphere<double,2>>(Sphere<double,2>(), Rect<double,1>(0,1));
 */
template <typename Shape>
class Extrusion : public virtual Convex<typename Shape::elem_t, Shape::N + 1> {
    public:
        /// The coordinate type of this Shape
        typedef typename Shape::elem_t T;
        /// The dimension of this Shape
        static constexpr size_t N = Shape::N + 1;
        
        /// Cross-section of this extrusion at the h = 0 plane.
        Shape base;
        /// Height range of this extrusion.
        Rect<T,1> height;
        
        
        /// Construct a new extrusion with unit height.
        Extrusion():height(0,1) {}
        
        /// Construct a new extrusion with cross section `base` and unit height.
        explicit Extrusion(const Shape& base):
            base(base),
            height(0,1) {}
        
        /// Construct an extrusion with cross section `base` and height range `height`.
        Extrusion(const Shape& base, const Rect<T,1>& height):
                base(base), 
                height(height) {}
        
        /** @brief Construct an extrusion with cross section `base`, and height 
         * ranging between `h0` and `h1`.
         */     
        Extrusion(const Shape& base, T h0, T h1):
            base(base),
            height(std::min(h0, h1), std::max(lo, h1)) {}
        
        /**
         * Extrusion-point intersection test.
         * 
         * @param p A point.
         * @return `true` if `p` is on or inside this extrusion; `false` otherwise.
         */
        bool contains(Vec<T,N> p) const {
            if (not height.contains(p[N-1])) return false;
            return base.contains(p.template resized<N-1>());
        }
        
        
        Vec<T,N> convex_support(Vec<T,N> d) const {
            typedef PointType<T,N-1> Pt;
            
            Vec<T,N-1> d_ = d.template resized<N-1>();
            
            if (d_.isZero()) {
                // support is normal to the endcap
                // pick an arbitrary point on the endcap
                d_[0] = 1;
            }
            
            Vec<T,N-1> p0 = base.convex_support(d_);
            
            // top or bottom face?
            return Vec<T,N>(p0, (d[N-1] > 0) ? height.hi : height.lo);
        }
        
        Rect<T,N> bounds() const {
            return shape.bounds() * height;
        }

        // todo: trace()


}; // class extrusion


/**
 * @brief Convenience function to extrude the shape `s` between heights `h0` and `h1` by
 * wrapping `s` in the `Extrusion` template.
 *
 * @related Extrusion
 */
template <typename Shape>
inline Extrusion<Shape> extrude(
        const    Shape& s,
        typename Shape::elem_t h0,
        typename Shape::elem_t h1)
{
    return Extrusion<Shape>(s, h0, h1);
}
    
} // namespace geom

#endif	/* EXTRUSION */

