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

#include "shapedetail/SeparatingAxis.h"

namespace geom {
    
    // TODO: Would be better if boxes had no translation?
    //       makes representation more unique. might cheapen certain ops.
 
// base class for OrientedRect. Not for direct use.
// We specialize the derived class because some dimensions don't support
// box-box intersection tests, and should not have methods for them.
template <typename T, int N>
class OrientedRectBase : virtual public Bounded<T,N> {
    
    public:
        
        /// Un-transformed extents.
        Rect<T,N> box;
        /// Transformation orienting `box`.
        AffineTransform<T,N> xf; // transforms object points to world
        
        /// Construct an axis-aligned box of size 1, with one corner at the origin.
        OrientedRectBase() {}
        
        /// Construct an axis-aligned box from the given Rect.
        OrientedRectBase(const Rect<T,N> &box):box(box) {}
        
        /// Construct an oriented box from the given Rect and object-to-world transformation.
        OrientedRectBase(const Rect<T,N> &box, const AffineTransform<T,N> &xf):
                box(box),
                xf(xf) {}
        
        
        /// Obtain an axis-aligned bounding box for this region.
        Rect<T,N> bounds() {
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
        
        
        /// Return the sphere circumscribing this OrientedRect.
        Sphere<T,N> sphereBounds() const {
            Vec<T,N> p0 = xf * box.min();
            Vec<T,N> p1 = xf * box.max();
            Vec<T,N> diag = (p1 - p0) / 2;
            return Sphere<T,N>(diag + p0, diag.mag());
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
        
};

/**
 * @ingroup shape
 * @brief A Rect with arbitrary position and orientation.
 */
#ifdef PARSING_DOXYGEN
template <typename T, index_t N>
#else
template <typename T, index_t N, typename Enable>
#endif
class OrientedRect : public OrientedRectBase<T,N> {
    
    private:
        typedef OrientedRectBase<T,N> base_t;
    
    public:
        
         /// Construct an axis-aligned box of size 1, with one corner at the origin.
        OrientedRect() {}
        
        /// Construct an axis-aligned box from the given Rect.
        OrientedRect(const Rect<T,N> &box):base_t(box) {}
        
        /// Construct an oriented box from the given Rect and object-to-world transformation.
        OrientedRect(const Rect<T,N> &box, const AffineTransform<T,N> &xf):
                base_t(box,xf) {}
        
        
        /**
         * Test whether this OrientedRect overlaps an axis-aligned box.
         * 
         * Not available for `N` > 3.
         * 
         * @param b1 Axis-aligned box to test.
         * @return `true` if and only if `this` overlaps with `b`; `false` otherwise. 
         */
        bool intersects(const Rect<T,N> &b1) {
            // here we apply the separating axis theorem.
            // all edge axes, face normals, and edge-edge normals must be tested.
            // the first two cases are one and the same; they are the principal
            // axes of each rect; amounting to 2N tests. The third case is
            // 9 axes for N=3 and 0 axes for N=2. Other dimensions may be
            // supported by extending detail::RectAxisHelper to higher N.
            const static index_t n_corners = 1 << N;
            Vec<T,N> b0_pts[n_corners];
            Vec<T,N> b1_pts[n_corners];
            Vec<T,N> b0_body_extreme[2] = { base_t::box.min(), base_t::box.max() };
            Vec<T,N> b1_body_extreme[2] = { b1.min(), b1.max()};
            
            // compute world-space points
            for (index_t c = 0; c < n_corners; c++) {
                for (index_t i = 0; i < N; i++) {
                    int min_or_max = (c & (1 << i)) != 0;
                    b0_pts[c][i] = b0_body_extreme[min_or_max][i];
                    b1_pts[c][i] = b1_body_extreme[min_or_max][i];
                }
                b0_pts[c] = base_t::xf * b0_pts[c];
            }
            
            // compare along world (i.e. b1's) axes
            // special case because no projection needed.
            for (index_t i = 0; i < N; i++) {
                T lo = std::numeric_limits<T>::max();
                T hi = std::numeric_limits<T>::lowest();
                for (index_t j = 0; j < n_corners; j++) {
                    lo = std::min(lo, b0_pts[j][i]);
                    hi = std::max(hi, b0_pts[j][i]);
                }
                // if no overlap on this axis, we've found a separating axis.
                if (lo >= b1.max()[i] or hi < b1.min()[i]) return false;
            }
            int m = 0;
            // now test all the remaining axes. For N=2, this is simply
            // the two principal axes of the oriented box. For N=3, we must also
            // test each axis which is the cross product of two principal axes,
            // for up to 12 additional axis tests (we have covered N already).
            for (detail::RectAxisHelper<T,N> helper(base_t::xf); not helper.done(); helper.next()) {
                m += 1;
                Vec<T,N> axis = helper.axis();
                if (axis == Vec<T,N>::zeros) continue; // redundant axis
                Rect<T,1> b0_range, b1_range;
                b0_range = b1_range = Rect<T,1>(std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest());
                for (index_t corner = 0; corner < n_corners; corner++) {
                    T b0_proj = axis.dot(b0_pts[corner]);
                    T b1_proj = axis.dot(b1_pts[corner]);
                    b0_range.extendTo(b0_proj);
                    b1_range.extendTo(b1_proj);
                }
                // no overlap on this axis; separating axis found.
                if (not b0_range.intersects(b1_range)) return false;
            }
            
            // overlap on all candidate axes.
            return true;
        }
        
        
        /**
         * Test whether this OrientedRect overlaps another.
         * 
         * Not available for `N` > 3.
         * @param b1 OritentedRect to test against.
         * @return `true` if and only if `this` overlaps with `b1`; `false` otherwise.
         */
        bool intersects(const OrientedRect<T,N> &b1) {
            // we will make ourselves axis-aligned; this is the same as applying
            // the inverse of our xf. We apply the inverse of xf to b1 too,
            // to preserve the relationhip between the us. From this, we
            // fallback to an ORect <-> Rect test.
            OrientedRect<T,N> b1_in_b0 = b1 / base_t::xf;
            return b1_in_b0.intersects(base_t::box);
        }

}; // OrientedRect


/**********************************************
 * OrientedRect case for when Separating Axis *
 * Test is not implemented; do not provide    *
 * intersection-test functions.               *
 **********************************************/


template <typename T, index_t N>
class OrientedRect<T,N, typename boost::enable_if_c<
                            not detail::RectAxisHelper<T,N>::SAT_SUPPORTED, 
                            void>::type > 
                       : public OrientedRectBase<T,N> {
    private:
        typedef OrientedRectBase<T,N> base_t;
    
    public:
        
        OrientedRect() {}
        
        OrientedRect(const Rect<T,N> &box):base_t(box) {}
        
        OrientedRect(const Rect<T,N> &box, const AffineTransform<T,N> &xf):
                base_t(box,xf) {}
};


/**********************************************
 * Operators                                  *
 **********************************************/

/** @addtogroup shape 
 *  @{
 */

/// Return a copy of `box` transformed by `xf`.
template <typename T, index_t N>
inline OrientedRect<T,N> operator*(const AffineTransform<T,N> &xf, OrientedRect<T,N> box) {
    box.xf *= xf;
    return box;
}


/// Return a copy of `box` transformed by the inverse of `xf`.
template <typename T, index_t N>
inline OrientedRect<T,N> operator/(OrientedRect<T,N> box, const AffineTransform<T,N> &xf) {
    box.xf /= xf;
    return box;
}
 

/// Apply `xf` to `b`.
template <typename T, index_t N>
inline OrientedRect<T,N>& operator*=(OrientedRect<T,N> &b, const AffineTransform<T,N> &xf) {
    b.xf *= xf;
    return b;
}

/// Apply the inverse of `xf` to `b`.
template <typename T, index_t N>
inline OrientedRect<T,N>& operator/=(OrientedRect<T,N> &b, const AffineTransform<T,N> &xf) {
    b.xf /= xf;
    return b;
}

/// @} // addtogroup linalg

} // namespace geom

#endif	/* ORIENTEDRECT_H */

