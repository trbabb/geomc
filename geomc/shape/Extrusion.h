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
 * @brief An N-dimensional extrusion of an arbitrary N-1 dimensional convex shape.
 * 
 * The first N-1 dimensions have cross-sections which are `ConvexShape`s
 * (e.g., `Rect`, `Sphere`, `Frustum`...). The shape is lofted "vertically" into the last dimension.
 *
 * Example:
 * 
 *     // initialize a cylinder of length and radius 1, centered at the origin:
 *     Extrusion<double, 3, Sphere> cylinder(Sphere<double,2>(), Rect<double,1>(0,1));
 */
template <typename T, index_t N, template <typename, index_t, typename... Args> ConvexShape>
class Extrusion : public virtual Convex<T,N> {
    public:
        /// Tranform orienting this extrusion.
        AffineTransform<T,N> xf;
        /// Cross-section of this extrusion at the h = 0 plane.
        ConvexShape<T,N-1,Args...> base;
        /// Height range of this extrusion.
        Rect<T,1> height;
        
        
        /// Construct a new extrusion, with unit length directly along the final ("height") axis.
        Extrusion():height(0,1) {}
        
        /// Construct an extrusion with cross section `base` and height range `height`.
        Extrusion(const ConvexShape<T,N-1,Args...>& base, const Rect<T,1>& height):
                base(base), 
                height(height) {}

        /// Construct an extrusion with cross section `base`, sheared axis `axis`, and height range `height`.
        Extrusion(const ConvexShape<T,N-1,Args...>& base, const Vec<T,N>& axis, const Rect<T,1>& height):
                base(base), 
                height(height) {
            // make it a pure shear:
            axis /= axis[N-1];
            // inverse of a shear is a shear in the opposite direction:
            Vec<T,N> inv = -axis;
            inv[N-1] = 1;
            for (index_t i = 0; i < N; ++i) {
                xf.mat[N-1][i] = axis[i];
                xf.inv[N-1][i] = inv[i];
            }
        }
        
        /**
         * Extrusion-point intersection test.
         * 
         * @param p A point.
         * @return `true` if `p` is on or inside this extrusion; `false` otherwise.
         */
        bool contains(Vec<T,N> p) const {
            p /= xf;
            if (not height.contains(p[N-1])) return false;
            return base.contains(p.template resized<N-1>());
        }
        
        
        Vec<T,N> convexSupport(Vec<T,N> d) const {
            typedef PointType<T,N-1> Pt;
            
            d = xf.applyInverseNormal(d);
            
            T hmax = height.max();
            T hmin = height.min();

            Vec<T,N-1> d_ = d.template resized<N-1>();

            if (d_.isZero()) {
                // support is normal to the endcap
                // pick an arbitrary point on the endcap
                d_[0] = 1;
            }
            
            Vec<T,N-1> p0 = base.convexSupport(d_);
            
            // top or bottom face?
            Vec<T,N> p1 = Vec<T,N>(p0, (d[N-1] > 0) ? hmax : hmin);
            
            return xf * p1;
        }

        // todo: trace()


}; // class extrusion
    
} // namespace geom

#endif	/* EXTRUSION */

