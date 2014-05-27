/* 
 * File:   SeparatingAxis.h
 * Author: tbabb
 *
 * Created on May 26, 2014, 11:27 PM
 */

#ifndef SEPARATINGAXIS_H
#define	SEPARATINGAXIS_H

/************************************************
 * Classes for iterating over axes to test for  *
 * separation when applying the separating axis *
 * theorem to test OrientedBox intersection.    *
 * The cases get far more complicated as N      *
 * increases. Only N=2 and N=3 are implemented. *
 ************************************************/

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
                    a[k] = xf.mat[j][k];
                }
                return a;
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
        
        // if arbitrary N are allowed, OrientedRect may become
        // a flat class, instead of the ugly enable_if<> business intended
        // to hide intersection test functions from higher dimensions.
        
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
                return this->_getXfAxis(i);
            } else { // we are assuming `not done()`
                // test all combinations of principal axes between
                // the two boxes. We assume one is axis-aligned.
                Vec<T,3> b0_axis;
                Vec<T,3> b1_axis = this->_getXfAxis(i % 3);
                b0_axis[i / 3] = 1;
                return b0_axis ^ b1_axis;
            }
        }
        
        inline bool done() {
            return base_t::i < 12; // 3 Xf axes, plus 3 * 3 axis combinations.
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
            return this->_getXfAxis(base_t::i);
        }
        
        inline bool done() {
            return base_t::i < 2;
        }
    };
    
} // namespace detail
    
} // namespace geom


#endif	/* SEPARATINGAXIS_H */

