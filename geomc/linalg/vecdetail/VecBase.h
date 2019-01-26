/*
 * VecBase.h
 *
 *  Created on: Oct 10, 2010
 *      Author: tbabb
 */

#ifndef VECBASE_H_
#define VECBASE_H_

#include <cmath>
#include <algorithm>

#if __cplusplus < 201103L
#include <tr1/functional>
#else
#include <functional>
#include <type_traits>
#endif

#include <boost/integer.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <geomc/Hash.h>
#include <geomc/linalg/LinalgTypes.h>

#ifdef GEOMC_VEC_CHECK_BOUNDS
#include <stdexcept>
#endif

#ifdef GEOMC_LINALG_USE_STREAMS
#include <iostream>
#endif

//TODO: specialize for dynamic dimension Vec<T,0>

//TODO: pare out neg(), add(), etc. which are holdovers from
//      javaland; also abs(), floor(), ceil(), etc. operators and
//      std:: functions are the preferred method. 
//      - make the std functions friends if possible, so that they show up 
//        with Vec in the docs?

//TODO: use c++11's decltype() and <auto> functionality to do proper widening/narrowing conversion
//      automatically, per http://stackoverflow.com/questions/9368432/c-arithmetic-operator-overloadingautomatic-widening
//      see also: boost::common_type<T,U>

//TODO: use begin() and end() pointers and pointer iteration instead of a counter?

namespace geom {
    
namespace detail {

    template <typename T, index_t N> class VecBase {
    public:

        VecBase() {
            std::fill(begin(), end(), 0);
        }

        VecBase(T a) {
            std::fill(begin(), end(), a);
        }

        VecBase(const T a[N]) {
            std::copy(a, a+N, begin());
        }
        
        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A const reference to the element at `idx`.
         */
        inline const T& get(index_t idx) const {
#ifdef GEOMC_VEC_CHECK_BOUNDS
                // this is considerably slower
                if (idx < 0 or idx >= N) {
                    throw std::out_of_range("vector index");
                }
#endif
            return v[idx];
        }
        
        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A reference to the element at `idx`.
         */
        inline T& get(index_t idx) {
#ifdef GEOMC_VEC_CHECK_BOUNDS
                if (idx < 0 or idx >= N) {
                    throw std::out_of_range("vector index");
                }
#endif
            return v[idx];
        }
        
        /// @return A read-only iterator pointing at the first element.
        inline const T* begin() const {
            return v;
        }
        
        /// @return A read-only iterator pointing just beyond the last element.
        inline const T* end() const {
            return v+N;
        }
        
        /// @return A writeable iterator pointing at the first element.
        inline T* begin() {
            return v;
        }
        
        /// @return A writeable iterator pointing just beyond the last element.
        inline T* end() {
            return v+N;
        }
        
    protected:
        T v[N];
    };
    
    template <typename T> class VecBase<T,2> {
    public:

        VecBase():x(0),y(0) {}

        VecBase(T a):x(a),y(a) {}

        VecBase(const T a[2]):x(a[0]),y(a[1]) {}
        
        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A const reference to the element at `idx`.
         */
        inline const T& get(index_t idx) const {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch(idx) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                throw std::out_of_range("vector index");
            }
#else
            return idx ? y : x;
#endif
        }
        
        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A reference to the element at `idx`.
         */
        inline T& get(index_t idx) {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch(idx) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                throw std::out_of_range("vector index");
            }
#else
                return idx ? y : x;
#endif
        }
        
        /// @return A read-only iterator pointing at the first element.
        inline const T* begin() const {
            return &x;
        }
        
        /// @return A read-only iterator pointing just beyond the last element.
        inline const T* end() const {
            return &y + 1;
        }
        
        /// @return A writeable iterator pointing at the first element.
        inline T* begin() {
            return &x;
        }
        
        /// @return A writeable iterator pointing just beyond the last element.
        inline T* end() {
            return &y + 1;
        }
        
        union {
            T x;
            T s;
            T row; // for matrix coordinates
        };
        
        union {
            T y;
            T t;
            T col;
        };
    };
    
    template <typename T> class VecBase<T,3> {
    public:

        VecBase():x(0),y(0),z(0) {}

        VecBase(T a):x(a),y(a),z(a) {}

        VecBase(const T a[3]):x(a[0]),y(a[1]),z(a[2]) {}

        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A const reference to the element at `idx`.
         */
        inline const T& get(index_t idx) const {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch (idx) {
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                default:
                    throw std::out_of_range("vector index"); //throw is very expensive, even if we don't use it. consider avoiding this.
            }
#else
            return (&x)[idx]; //faster than a 4-branch case, whether the case has a throw or not (gcc)
#endif
        }
        
        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A reference to the element at `idx`.
         */
        inline T& get(index_t idx) {
#ifdef GEOMC_CHECK_BOUNDS
            switch (idx) {
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                default:
                    throw std::out_of_range("vector index"); //very slow; whether used or not.
            }
#else
            return (&x)[idx];
#endif
        }
        
        /// @return A read-only iterator pointing at the first element.
        inline const T* begin() const {
            return &x;
        }
        
        /// @return A read-only iterator pointing just beyond the last element.
        inline const T* end() const {
            return &z + 1;
        }
        
        /// @return A writeable iterator pointing at the first element.
        inline T* begin() {
            return &x;
        }
        
        /// @return A writeable iterator pointing just beyond the last element.
        inline T* end() {
            return &z + 1;
        }
        
        union {
            T x;
            T s;
            T r;
        };
        union {
            T y;
            T t;
            T g;
        };
        union {
            T z;
            T u;
            T b;
        };
    };
    
    template <typename T> class VecBase<T,4> {
    public:

        VecBase():x(0),y(0),z(0),w(0) {}

        VecBase(T a):x(a),y(a),z(a),w(a) {}

        VecBase(const T a[4]):x(a[0]),y(a[1]),z(a[2]),w(a[3]) {}

        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A const reference to the element at `idx`.
         */
        inline const T& get(index_t idx) const {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch (idx) {
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                case 3:
                    return w;
                default:
                    throw std::out_of_range("vector index"); //very slow; whether used or not.
            }
#else
            return (&x)[idx];
#endif
        }
        
        /**
         * Get the element at index `idx`.
         * @param idx Index of element.
         * @return A reference to the element at `idx`.
         */
        inline T& get(index_t idx) {
#ifdef GEOMC_VEC_CHECK_BOUNDS
            switch (idx) {
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                case 3:
                    return w;
                default:
                    throw std::out_of_range("vector index");
            }
#else
            return (&x)[idx];
#endif
        }
        
        /// @return A read-only iterator pointing at the first element.
        inline const T* begin() const {
            return &x;
        }
        
        /// @return A read-only iterator pointing just beyond the last element.
        inline const T* end() const {
            return &w + 1;
        }
        
        /// @return A writeable iterator pointing at the first element.
        inline T* begin() {
            return &x;
        }
        
        /// @return A writeable iterator pointing just beyond the last element.
        inline T* end() {
            return &w + 1;
        }
        
        union {
            T x;
            T s;
            T r;
        };
        union {
            T y;
            T t;
            T g;
        };
        union {
            T z;
            T u;
            T b;
        };
        union {
            T w;
            T v;
            T a;
        };
    };
    
    // end VecBase classes
    
    
    template <typename T, index_t N, typename VType> class VecCommon : public VecBase<T,N> {
        public:
        
        /// Self type. I.e., `Vec<T,N>` if a vector, `Quat<T>` if a quaternion.
        typedef VType self_t;
        /// Element type.
        typedef T elem_t;
        /// Vector dimension.
        static const index_t DIM = N;
        
        static const self_t ones;
        static const self_t zeros;

        /*******************************
         * Structors                   *
         *******************************/

        VecCommon():VecBase<T,N>() {}

        VecCommon(T a):VecBase<T,N>(a) {}

        VecCommon(const T a[N]):VecBase<T,N>(a) {}

        /*******************************
         * Operators                   *
         *******************************/

        /**
         * Vector element access.
         * 
         * @param idx Index of element to retrieve.
         * @return A read-only reference to the element at index `idx`.
         */
        inline const T& operator[](index_t idx) const {
            return this->get(idx);
        }
        
        /**
         * Vector element access.
         * 
         * @param idx Index of element to retrieve.
         * @return A reference to the element at index `idx`.
         */
        inline T& operator[](index_t idx) {
            return this->get(idx);
        }

        /**
         * Element typecast.
         * 
         * @return A new vector with all elements cast to type `U`.
         */
        template <typename U> 
        explicit operator Vec<U,N>() const {
            Vec<U,N> r;
            for (index_t i = 0; i < N; i++) {
                r[i] = (U)(this->get(i));
            }
            return r;
        }

        /// @return `true` if all corresponding elements of `this` and `vv` are equal, `false` otherwise.
        inline bool operator==(const self_t &vv) const {
            for (index_t i = 0; i < N; i++) {
                if (vv.get(i) != this->get(i)) return false;
            }
            return true;
        }
        
        /// @return `true` if any corresponding elements of `this` and `vv` are unequal, `false` otherwise.
        inline bool operator!=(const self_t &vv) const {
            for (index_t i = 0; i < N; i++) {
                if (vv.get(i) != this->get(i)) return true;
            }
            return false;
        }

        /// Element-wise addition.
        inline self_t operator+(const self_t &v) const {
            return this->add(v);
        }

        /// Element-wise addition and assignment.
        self_t& operator+=(const self_t &vv) {
            for (index_t i = 0; i < N; i++) {
                this->get(i) += vv[i];
            }
            return *(static_cast<self_t*>(this));
        }
        
        /// Element-wise subtraction.
        inline self_t operator-(const self_t &v) const {
            return this->sub(v);
        }
        
        /// Subtraction and assignment.
        self_t& operator-=(const self_t &vv) {
            for (index_t i = 0; i < N; i++) {
                this->get(i) -= vv[i];
            }
            return *(static_cast<self_t*>(this));
        }

        /// Scalar multiplication and assignment.
        self_t& operator*=(T s) {
            for (index_t i = 0; i < N; i++) {
                this->get(i) *= s;
            }
            return *(static_cast<self_t*>(this));
        }

        /// Scalar division and assignment.
        self_t& operator/=(T s) {
            for (index_t i = 0; i < N; i++) {
                this->get(i) /= s;
            }
            return *(static_cast<self_t*>(this));
        }

        /// Element-wise multiplication and assignment.
        self_t& operator*=(const self_t &vv) {
            for (index_t i = 0; i < N; i++) {
                this->get(i) *= vv[i];
            }
            return *(static_cast<self_t*>(this));
        }

        /** @return A copy of this vector with all elements negated (i.e. a 
         * vector pointing in the opposite direction).
         */
        inline self_t operator-() const {
            return this->neg();
        }

        /*******************************
         * Methods                     *
         *******************************/
        
        /**
         * @param v Another vector.
         * @return A new vector `x` such that `x[i] = this[i] + v[i]`.
         */
        inline self_t add(const self_t &v) const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = this->get(i) + v.get(i);
            }
            return r;
        }
        
        /**
         * @param v Another vector.
         * @return A new vector `x` such that `x[i] = this[i] - v[i]`.
         */
        inline self_t sub(const self_t &v) const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = this->get(i) - v[i];
            }
            return r;
        }
        
        /**
         * @param v Another vector.
         * @return A new vector `x` such that `x[i] = this[i] * v[i]`.
         */
        inline self_t scale(const self_t &v) const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = this->get(i) * v[i];
            }
            return r;
        }
        
        /**
         * @param a A constant scale factor.
         * @return A new vector `x` such that `x[i] = this[i] * a`.
         */
        inline self_t scale(T a) const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = this->get(i) * a;
            }
            return r;
        }
        
        /**
         * @return A new vector `x` such that `x[i] = -this[i]`.
         */
        inline self_t neg() const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = -this->get(i);
            }
            return r;
        }
 
        
        /**
         * @return A copy of this vector with unit length.
         */
        inline self_t unit() const {
            T n = mag();
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = this->get(i) / n;
            }
            return r;
        }
        
        
        /**
         * @param v Another vector.
         * @return The dot product of `this` with `v`.
         */
        inline T dot(const self_t &v) const {
            T sum = 0;
            for (index_t i = 0; i < N; i++) {
                sum += this->get(i) * v[i];
            }
            return sum;
        }
        
        /**
         * @return The magnitude (geometric length) of this vector.
         */
        inline T mag() const {
            return std::sqrt(this->mag2());
        }
        
        /**
         * @return The square of the magnitude of this vector.
         */
        inline T mag2() const {
            T sum = 0;
            for (index_t i = 0; i < N; i++) {
                T cur = this->get(i);
                sum += cur * cur;
            }
            return sum;
        }

        /**
         * @param pt Another point.
         * @return The distance between `this` and `pt`.
         */
        T dist(const self_t &pt) const {
            return (this->sub(pt)).mag();
        }

        /**
         * @param pt Another point.
         * @return The square of the distance between `this` and `pt`.
         */
        T dist2(const self_t &pt) const {
            return (this->sub(pt)).mag2();
        }

        /**
         * @param axis Axis of reflection.
         * @return A copy of this vector reflected across the given axis.
         */
        self_t reflect(self_t axis) const {
            axis = axis.unit();
            return this->neg() + axis * 2 * this->dot(axis);
        }

        /**
         * Treat `this` as a velocity vector or incident ray; this function returns
         * the velocity reflected off of a surface with normal `normal`. 
         * Equivalent to `-reflect(normal)`.
         * 
         * @param normal Normal of surface to "bounce" on.
         * @return The "bounced" direction vector.
         */
        self_t bounce(const self_t &normal) const {
            return -reflect(normal);
        }
        
        /**
         * @param v A direction vector.
         * @return A vector in direction `v` with magnitude equal to the component
         * of `this` aligned with `v`.
         */
        self_t projectOn(const self_t &v) const {
            T m2 = v.mag2();
            return m2 != 0 ? (v * this->dot(v) / m2) : v;
        }

        /**
         * Linear interpolation. A mix parameter of 0 evaluates to `this`, while 
         * 1 is `v`.
         * @param v Another vector.
         * @param mix A mixing factor between 0 and 1.
         * @return A linear mixing of `this` with `v`.
         */
        self_t mix(const self_t &v, T mix) const {
            T xim = 1 - mix;
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++) {
                r[i] = xim*this->get(i) + mix*v[i];
            }
            return r;
        }

        /**
         * @param v Another vector.
         * @return Angle in radians between `this` and `v`, between 0 and `pi`.
         */
        T angleTo(const self_t &v) const {
            return std::acos(v.unit().dot(unit()));
        }
        
        ///////////////////////
        // Clamping/Rounding //
        ///////////////////////
        
        /**
         * Element-wise absolute value.
         * @return A new vector `x` such that `x[i] = abs(this[i])`.
         */
        self_t abs() const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = std::abs(this->get(i));
            }
            return r;
        }
        
        //todo: handle integer cases
        /**
         * Element-wise floor function.
         * @return A new vector `x` such that `x[i] = floor(this[i])`.
         */
        self_t floor() const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = std::floor(this->get(i));
            }
            return r;
        }
        
        //todo: handle integer cases
        /**
         * Element-wise ceiling function.
         * @return A new vector `x` such that `x[i] = ceil(this[i])`.
         */
        self_t ceil() const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = std::ceil(this->get(i));
            }
            return r;
        }
        
        /**
         * Element-wise minimum.
         * @param v Another vector.
         * @return A new vector `x` such that `x[i] = min(this[i], v[i])`.
         */
        self_t min(const self_t &v) const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = std::min(this->get(i), v[i]);
            }
            return r;
        }
        
        /**
         * Element-wise maximum.
         * @param v Another vector.
         * @return A new vector `x` such that `x[i] = max(this[i], v[i])`.
         */
        self_t max(const self_t &v) const {
            self_t r;
            for (index_t i = 0; i < N; i++) {
                r[i] = std::max(this->get(i), v[i]);
            }
            return r;
        }
        
        /**
         * Element-wise clamp.
         * @param lo Element-wise lower extremes.
         * @param hi Element-wise upper extremes.
         * @return A new vector such that each element `x[i]` is clamped
         * between `lo[i]` and `hi[i]`.
         */
         self_t clamp(const self_t &lo, const self_t &hi) const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++) {
                r[i]  = std::min(std::max(this->get(i), lo[i]), hi[i]);
            }
            return r;
        }
        
        /**
         * @return A new vector with each element rounded to the nearest integer.
         */
        self_t round() const {
            Vec<T,N> r;
            for (index_t i = 0; i < N; i++) {
                r[i] = std::floor(this->get(i) + 0.5);
            }
            return r;
        }
        
        /**
         * Resized copy of a vector.
         * @tparam M Dimension of new vector.
         * @return A new vector with size `M`. If `M` is largner than `N`, the
         * new elements will be set to zero.
         */
        template <index_t M> inline Vec<T,M> resized() const {
            Vec<T,M> out;
            const T *p = this->begin();
            std::copy(p, p + std::min(M,N), out.begin());
            return out;
        }

        /**
         * @return `true` if all elements are zero, `false` otherwise.
         */
        bool isZero() const {
            for (index_t i = 0; i < N; i++) {
                if (this->get(i) != (T)0) return false;
            }
            return true;
        }

        /**
         * @return The number of elements in this vector. Always equal to `N`. 
         */
        inline index_t getSize() const {
            return N;
        }
        
        /**
         * @return A set of pseudo-random bits deterministically related to the
         * elements of this vector.
         */
        inline size_t hashcode() const {
            return general_hash(this, sizeof(Vec<T,N>));
        }

    }; /* end VecCommon class */

    template <typename T, index_t N, typename VType> const VType VecCommon< T,N,VType >::ones  = VType((T)1);
    template <typename T, index_t N, typename VType> const VType VecCommon< T,N,VType >::zeros = VType((T)0);

    // used by fromRGB()
    // for integer vectors, you don't want to divide by 255, essentially
    // these perform that logic switch:
    
    template <typename T, typename Enable=void>
    struct _RGBChannelConversion {
        static const double scale;
    };
    
    // technically c++ forbids use of floats in constant expressions, for some weird reason.
    // gcc does not complain unless compiling under c++11, ironically enough where
    // the requirements for constant expressions are *supposed* to be relaxed.
    template <typename T, typename Enable>
    const double _RGBChannelConversion<T,Enable>::scale = 1.0 / 255.0;
    
    template <typename T>
    struct _RGBChannelConversion<T, typename boost::enable_if<boost::is_integral<T>,void>::type> {
        static const int scale = 1;
    };
    
    template <typename T, typename Enable=void>
    struct IsVector {
        static const bool value = false;
    };
    
    template <typename T>
    struct IsVector<T, typename boost::enable_if<
                            boost::is_base_of<
                                detail::VecCommon<typename T::elem_t, T::DIM, T>, 
                                T>,
                        void>::type
                     > {
        const static bool value = true;
    };
    
} /* end namespace detail */

} /* end namepsace geom */


/*********************************
 * Hashing                       *
 *********************************/

//allows Vec<T,N>s to work with std::tr1 hash containers
namespace std { 

#if __cplusplus < 201103L
namespace tr1 {
#endif
    
   template <typename T, index_t N> struct hash< geom::Vec<T,N> > : public unary_function<geom::Vec<T,N>, size_t> {
       size_t operator()(const geom::Vec<T,N>& v) const {
           return v.hashcode();
       }
   };

#if __cplusplus < 201103L
}
#endif

}

#endif /* VECBASE_H_ */
