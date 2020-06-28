/* 
 * File:   SeparatingAxis.h
 * Author: tbabb
 *
 * Created on May 26, 2014, 11:27 PM
 */

#ifndef SEPARATINGAXIS_H
#define	SEPARATINGAXIS_H


/************************************************************
 * For small dimensions, we use the separating axis theorem *
 * to perform box-box intersection tests. The number of     *
 * axes we must test and the strategy for choosing them     *
 * changes based on the dimension; these classes are        *
 * specialized for those dimensions. The cases get far more *
 * complicated as N increases; only N=2 and N=3 are         *
 * implemented.                                             *
 ************************************************************/


namespace geom {
    
namespace detail {
    
    
    /***************************************
     * Iterator base class                 *
     ***************************************/
    
    template <typename T, index_t N>
    class RectAxisHelperBase {
        
        public:
        
            const AffineTransform<T,N> &xf;
            index_t i;

            RectAxisHelperBase(const AffineTransform<T,N> &xf) : xf(xf), i(0) {}

            inline void next() {
                i += 1;
            }
        
        protected:
            
            Vec<T,N> _getXfAxis(index_t j) {
                Vec<T,N> a;
                for (index_t k = 0; k < N; k++) {
                    a[k] = xf.mat[k][j];
                }
                return a;
            }
        
            
            Vec<T,N> _getNormal(index_t j) {
                Vec<T,N> n;
                n[j] = 1;
                return xf.applyNormal(n);
            }
        
    };
    
    
    /***************************************
     * General N-dimensional axis tester.  *
     * Not supported for N != {2,3}.       *
     ***************************************/
    
    template <typename T, index_t N>
    class RectAxisHelper : public RectAxisHelperBase<T,N> {
        public: 
        static const bool SAT_SUPPORTED = false;
        
        // add support for arbitrary dimensions here.
        // (GJK may be best for the general case, however)
        
    };
    
    
    /********** N=3 Axis Iterator **********/
    
    template <typename T>
    class RectAxisHelper<T,3> : public RectAxisHelperBase<T,3> {
        private:
            
            typedef RectAxisHelperBase<T,3> base_t;
            
        public:
        
        static const bool SAT_SUPPORTED = true;
        
        RectAxisHelper(const AffineTransform<T,3> &xf) : base_t(xf) {}
        
        // we don't include the canonical axes, because we handle
        // them separately for speed. A more user-facing implementation
        // would iterate over all 15 axes. 
        
        Vec<T,3> axis() {
            index_t i = base_t::i;
            if (i < 3) {
                // test each of the principal axes in the oriented box
                return this->_getNormal(i);
            } else { // we are assuming `not done()`
                // test all combinations of principal axes between
                // the two boxes. We assume one is axis-aligned.
                Vec<T,3> b0_axis;
                Vec<T,3> b1_axis = this->_getXfAxis(i % 3);
                b0_axis[i / 3 - 1] = 1;
                return b0_axis ^ b1_axis;
            }
        }
        
        inline bool done() {
            return base_t::i >= 12; // 3 Xf axes, plus 3 * 3 axis combinations.
        }
    };
    
    
    /********** N=2 Axis Iterator **********/
    
    template <typename T>
    class RectAxisHelper<T,2> : public RectAxisHelperBase<T,2> {
        private:
            typedef RectAxisHelperBase<T,2> base_t;
        public:
        
        static const bool SAT_SUPPORTED = true;
        
        RectAxisHelper(const AffineTransform<T,2> &xf) : base_t(xf) {}
        
        // canonical axes already handled.
        
        Vec<T,2> axis() {
            return this->_getNormal(base_t::i);
        }
        
        inline bool done() {
            return base_t::i >= 2;
        }
    };
    
    
    /***********************************************************
     * Intersector class. Use separating axis theorem if       *
     * supported; else gjk.                                    *
     ***********************************************************/
    
    template <typename T, index_t N, typename Enable=void>
    class RectIntersector {
        
        static inline bool intersect(const OrientedRect<T,N> &b0, const Rect<T,N> &b1) {
            return gjk_intersect(b0, b1);
        }
    };
    
    template <typename T, index_t N>
    class RectIntersector<T, N, typename boost::enable_if_c<RectAxisHelper<T,N>::SAT_SUPPORTED>::type> {
        public: 
            
        static bool intersect(const OrientedRect<T,N> &b0, const Rect<T,N> &b1) {
            // here we apply the separating axis theorem.
            // all edge axes, face normals, and edge-edge normals must be tested.
            // the first two cases are one and the same; they are the principal
            // axes of each rect; amounting to 2N tests. The third case is
            // 9 axes for N=3 and 0 axes for N=2. Other dimensions may be
            // supported by extending detail::RectAxisHelper to higher N.
            const static index_t n_corners = 1 << N;
            Vec<T,N> b0_pts[n_corners];
            Vec<T,N> b1_pts[n_corners];
            Vec<T,N> b0_body_extreme[2] = { b0.box.lo, b0.box.hi };
            Vec<T,N> b1_body_extreme[2] = {     b1.lo,     b1.hi };
            
            // compute world-space points
            for (index_t c = 0; c < n_corners; c++) {
                for (index_t i = 0; i < N; i++) {
                    int min_or_max = (c & (1 << i)) != 0;
                    b0_pts[c][i] = b0_body_extreme[min_or_max][i];
                    b1_pts[c][i] = b1_body_extreme[min_or_max][i];
                }
                b0_pts[c] = b0.xf * b0_pts[c];
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
                if (lo >= b1.hi[i] or hi < b1.lo[i]) return false;
            }
            
            // now test all the remaining axes. For N=2, this is simply
            // the two principal axes of the oriented box. For N=3, we must also
            // test each axis which is the cross product of two principal axes,
            // for up to 12 additional axis tests (we have covered N already).
            for (detail::RectAxisHelper<T,N> helper(b0.xf); not helper.done(); helper.next()) {
                Vec<T,N> axis = helper.axis();
                if (axis == Vec<T,N>::zeros) continue; // redundant axis
                Rect<T,1> b0_range, b1_range;
                b0_range = b1_range = Rect<T,1>(std::numeric_limits<T>::max(), 
                                                std::numeric_limits<T>::lowest());
                for (index_t corner = 0; corner < n_corners; corner++) {
                    T b0_proj = axis.dot(b0_pts[corner]);
                    T b1_proj = axis.dot(b1_pts[corner]);
                    b0_range |= b0_proj;
                    b1_range |= b1_proj;
                }
                // no overlap on this axis; separating axis found.
                if (not b0_range.intersects(b1_range)) return false;
            }
            
            // overlap on all candidate axes.
            return true;
        }
        
    };
    
    
} // namespace detail
    
} // namespace geom


#endif	/* SEPARATINGAXIS_H */

