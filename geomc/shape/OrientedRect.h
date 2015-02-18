/* 
 * File:   OrientedRect.h
 * Author: tbabb
 *
 * Created on May 26, 2014, 6:58 PM
 */

#ifndef ORIENTEDRECT_H
#define	ORIENTEDRECT_H

#include <boost/utility/enable_if.hpp>
#include <geomc/linalg/AffineTransform.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/shapedetail/SeparatingAxis.h>

#include "Sphere.h"

namespace geom {
  
/**
 * @ingroup shape
 * @brief A Rect with arbitrary position and orientation.
 */
template <typename T, index_t N>
class OrientedRect : virtual public Bounded<T,N>, virtual public Convex<T,N> {
    
    public:
        
        /// Un-transformed extents.
        Rect<T,N> box;
        /// Transformation orienting `box`.
        AffineTransform<T,N> xf; // transforms object points to world
        
        /// Construct an axis-aligned box of size 1, with one corner at the origin.
        OrientedRect() {}
        
        /// Construct an axis-aligned box from the given Rect.
        OrientedRect(const Rect<T,N> &box):box(box) {}
        
        /// Construct an oriented box from the given Rect and object-to-world transformation.
        OrientedRect(const Rect<T,N> &box, const AffineTransform<T,N> &xf):
                box(box),
                xf(xf) {}
        
        
        /// Obtain an axis-aligned bounding box for this region.
        Rect<T,N> bounds() {
            // faster than using convexSupport, somewhat surprisingly.
            Vec<T,N> pts[2] = {box.min(), box.max()};
            Vec<T,N> lo = xf * box.min();
            Vec<T,N> hi = lo;
            
            // for each corner
            for (index_t i = 1; i < (1 << N); i++) {
                Vec<T,N> p;
                // for each coordinate
                for (index_t j = 0; j < N; j++) {
                    p[j] = pts[(i & (1 << j)) != 0][j];
                }
                p = xf * p;
                lo = std::min(lo, p);
                hi = std::max(hi, p);
            }
            
            return Rect<T,N>(lo,hi);
        }
        
        
        Vec<T,N> convexSupport(Vec<T,N> d) const {
            Vec<T,N> d_body = d / xf;
            Vec<T,N> o;
            for (index_t i = 0; i < N; i++) {
                o[i] = ((d_body[i] < 0) ? box.min() : box.max())[i];
            }
            return xf * o;
        }
        
        
        /**
         * @brief Fill an array with the vertecies of this OrientedBox.
         * 
         * @param out Array to receive 2<sup>N</sup> corner vertecies.
         */
        void getCorners(Vec<T,N> out[1 << N]) {
            const static index_t n_corners = 1 << N;
            Vec<T,N> extreme[2] = { box.min(), box.max() };
            for (index_t c = 0; c < n_corners; c++) {
                Vec<T,N> pt;
                for (index_t i = 0; i < N; i++) {
                    pt[i] = extreme[(c & (1 << i)) != 0][i];
                }
                out[c] = xf * pt;
            }
        }
        
        
        /// Returns true if and only if `p` is inside or on the surface of this OrientedRect.
        inline bool contains(Vec<T,N> p) {
            p = p / xf; // world coord to object coord
            return box.contains(p);
        }
        
        /**
         * Test whether this OrientedRect overlaps an axis-aligned Rect.
         * 
         * @param b1 OritentedRect to test against.
         * @return `true` if and only if `this` overlaps with `b1`; `false` otherwise.
         */
        inline bool intersects(Rect<T,N> r) {
            return detail::RectIntersector<T,N>::intersect(*this, r);
        }
        
        /**
         * Test whether this OrientedRect overlaps another.
         * 
         * @param b1 OritentedRect to test against.
         * @return `true` if and only if `this` overlaps with `b1`; `false` otherwise.
         */
        inline bool intersects(const OrientedRect<T,N> &b1) {
            // we will make ourselves axis-aligned; this is the same as applying
            // the inverse of our xf. We apply the inverse of xf to b1 too,
            // to preserve the relationhip between the us. From this, we
            // fallback to an ORect <-> Rect test.
            OrientedRect<T,N> b1_in_b0 = b1 / xf;
            
            return b1_in_b0.intersects(box);
        }
        
        
        /**
         * Box-ray intersection test.
         * 
         * @param r A ray
         * @param sides Whether to hit-test the front-facing or back-facing surfaces.
         * @return A ray Hit describing whether and where the ray intersects this OrientedRect,
         * as well as the normal, side hit, and ray parameter.
         */
        Hit<T,N> trace(const Ray<T,N> &r, HitSide side) {
            Ray<T,N> rbody = r / xf;
            Hit<T,N> h = box.trace(rbody, side);
            if (h.hit) {
                h.p = xf * h.p;
                h.n = xf.applyNormal(h.n);
            }
            return h;
        }
        
}; // class OrientedRect


/**********************************************
 * Operators                                  *
 **********************************************/


/** Return a copy of `box` transformed by `xf`.
 * @related OrientedRect
 */
template <typename T, index_t N>
inline OrientedRect<T,N> operator*(const AffineTransform<T,N> &xf, OrientedRect<T,N> box) {
    box.xf *= xf;
    return box;
}


/** Return a copy of `box` transformed by the inverse of `xf`.
 * @related OrientedRect
 */
template <typename T, index_t N>
inline OrientedRect<T,N> operator/(OrientedRect<T,N> box, const AffineTransform<T,N> &xf) {
    box.xf /= xf;
    return box;
}
 

/** Apply `xf` to `b`.
 * @related OrientedRect
 */
template <typename T, index_t N>
inline OrientedRect<T,N>& operator*=(OrientedRect<T,N> &b, const AffineTransform<T,N> &xf) {
    b.xf *= xf;
    return b;
}

/** Apply the inverse of `xf` to `b`.
 * @related OrientedRect
 */
template <typename T, index_t N>
inline OrientedRect<T,N>& operator/=(OrientedRect<T,N> &b, const AffineTransform<T,N> &xf) {
    b.xf /= xf;
    return b;
}


} // namespace geom

#endif	/* ORIENTEDRECT_H */

