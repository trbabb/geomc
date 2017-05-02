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
template <typename T, index_t N>
class Frustum : public virtual Convex<T,N> {
    public:
        /// Tranform orienting this frustum.
        AffineTransform<T,N> xf;
        /// Extent of this frustum at the `h = 1` plane.
        Rect<T,N-1> base;
        /// Height range spanned by this frustum.
        Rect<T,1> height;
        
        
        /// Construct a new frustum, spanning the unit rectange and heights between 1 and 2.
        Frustum(): base(typename Rect<T,N-1>::point_t(-1), 
                        typename Rect<T,N-1>::point_t( 1)), 
                   height(1,2) {}
        
        /// Construct a new frustum spanning the rectangle given by `base` and the height range between `h0` and `h1`.
        Frustum(Rect<T,N-1> base, T h0, T h1):
                    base(base), 
                    height(std::min(h0,h1), std::max(h0,h1)) {}
        
        /**
         * Frustum-point intersection test.
         * 
         * @param p A point.
         * @return `true` if `p` is on or inside this frustum; `false` otherwise.
         */
        bool contains(Vec<T,N> p) const {
            p /= xf;
            if (not height.contains(p[N-1])) return false;
            p /= p[N-1];
            return base.contains(p.template resized<N-1>());
        }
        
        
        Vec<T,N> convexSupport(Vec<T,N> d) const {
            typedef PointType<T,N-1> Pt;
            
            d = xf.applyInverseNormal(d);
            
            T hmax = height.max();
            T hmin = height.min();
            T sign = (hmax < 0) ? -1 : 1;
            
            // find the base corner
            Vec<T,N> p;
            for (index_t i = 0; i < N-1; i++) {
                if (sign * d[i] > 0) p[i] = Pt::iterator(base.max())[i];
                else                 p[i] = Pt::iterator(base.min())[i];
            }
            p[N-1] = 1;
            
            // top or bottom corner?
            p *= (p.dot(d) > 0) ? hmax : ((hmax > 0 and hmin < 0) ? 0 : hmin);
            
            return xf * p;
        }
        
        /**
         * Fill `p` with the positions of the 2<sup>N</sup> corners of
         * this frustum.
         * @param p Array to receive 2<sup>N</sup> corner positions.
         */
        void getCorners(Vec<T,N> p[1 << N]) const {
            typedef PointType<T,N-1> Pt;
            typename Pt::point_t extreme[2] = { base.min(), base.max() };
            
            const T hmax = height.max();
            T hmin = height.min();
            hmin = (hmax > 0 and hmin < 0) ? 0 : hmin;
            
            for (index_t k = 0; k < 2; k++) {
                // for lower and upper extents
                T h = k == 0 ? hmin : hmax;
                for (index_t c = 0; c < (1 << (N-1)); c++) {
                    // for each corner of the base rect
                    Vec<T,N> pt;
                    for (index_t i = 0; i < N-1; i++) {
                        // for each coordinate of the corner point
                        pt[i] = Pt::iterator(extreme[(c & (1 << i)) != 0])[i];
                    }
                    pt[N-1] = 1;
                    pt *= h;
                    pt  = xf * pt;
                    p[c + ((k > 0) ? (1 <<(N-1)) : 0)] = pt;
                }
            }
        }

}; // class frustum
    
} // namespace geom

#endif	/* FRUSTUM_H */

