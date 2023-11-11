#pragma once

/*
 * VecBase.h
 *
 *  Created on: Oct 10, 2010
 *      Author: tbabb
 */

#include <cmath>
#include <algorithm>
#include <functional>
#include <type_traits>

#include <geomc/linalg/LinalgTypes.h>

//TODO: specialize for dynamic dimension Vec<T,0>

//TODO: pare out neg(), add(), etc. which are holdovers from
//      javaland; also abs(), floor(), ceil(), etc. operators and
//      std:: functions are the preferred method. 
//      - make the std functions friends if possible, so that they show up 
//        with Vec in the docs?

//TODO: use c++11's decltype() and <auto> functionality to do proper widening/narrowing
//      conversion automatically, per 
//         http://stackoverflow.com/questions/9368432/
//            c-arithmetic-operator-overloadingautomatic-widening
//      see also: std::common_type<T,U>

//TODO: use begin() and end() pointers and pointer iteration instead of a counter?

namespace geom {
    
namespace detail {

template <typename T> 
inline bool is_real(T d){
    return !(std::isinf(d) || std::isnan(d));
}

template <typename T, index_t N>
class VecBase {
protected:

    constexpr VecBase() {
        for (index_t i = 0; i < N; ++i) {
            v[i] = 0;
        }
        // not constexpr until c++20:
        // std::fill(begin(), end(), 0);
    }

    constexpr VecBase(T a) {
        for (index_t i = 0; i < N; ++i) {
            v[i] = a;
        }
        // not constexpr until c++20:
        // std::fill(begin(), end(), a);
    }

    constexpr VecBase(const T a[N]) {
        for (index_t i = 0; i < N; ++i) {
            v[i] = a[i];
        }
        // not constexpr until c++20:
        // std::copy(a, a+N, begin());
    }

public:
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A const reference to the element at `idx`.
     */
    inline const T& get(index_t idx) const {
        return v[idx];
    }
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A reference to the element at `idx`.
     */
    inline T& get(index_t idx) {
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

template <typename T>
class VecBase<T,2> {
protected:

    constexpr VecBase():x(0),y(0) {}

    constexpr VecBase(T a):x(a),y(a) {}

    constexpr VecBase(const T a[2]):x(a[0]),y(a[1]) {}

public:
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A const reference to the element at `idx`.
     */
    inline const T& get(index_t idx) const {
        return idx ? y : x;
    }
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A reference to the element at `idx`.
     */
    inline T& get(index_t idx) {
        return idx ? y : x;
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

template <typename T>
class VecBase<T,3> {
protected:

    constexpr VecBase():x(0),y(0),z(0) {}

    constexpr VecBase(T a):x(a),y(a),z(a) {}

    constexpr VecBase(const T a[3]):x(a[0]),y(a[1]),z(a[2]) {}

public:
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A const reference to the element at `idx`.
     */
    inline const T& get(index_t idx) const {
        return (&x)[idx];
    }
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A reference to the element at `idx`.
     */
    inline T& get(index_t idx) {
        return (&x)[idx];
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

template <typename T>
class VecBase<T,4> {
protected:
    
    constexpr VecBase():x(0),y(0),z(0),w(0) {}

    constexpr VecBase(T a):x(a),y(a),z(a),w(a) {}

    constexpr VecBase(const T a[4]):x(a[0]),y(a[1]),z(a[2]),w(a[3]) {}

public:
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A const reference to the element at `idx`.
     */
    inline const T& get(index_t idx) const {
        return (&x)[idx];
    }
    
    /**
     * Get the element at index `idx`.
     * @param idx Index of element.
     * @return A reference to the element at `idx`.
     */
    inline T& get(index_t idx) {
        return (&x)[idx];
    }
    
    /// A read-only iterator pointing at the first element.
    inline const T* begin() const {
        return &x;
    }
    
    /// A read-only iterator pointing just beyond the last element.
    inline const T* end() const {
        return &w + 1;
    }
    
    /// A writeable iterator pointing at the first element.
    inline T* begin() {
        return &x;
    }
    
    /// A writeable iterator pointing just beyond the last element.
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


/// Common base for all Vec-derived classes. Do not instantiate directly.
template <typename T, index_t N, typename VType>
class VecCommon : public VecBase<T,N> {
public:
    
    static_assert(N > 0, "Vector dimension must be larger than zero and cannot be dynamic");
    
    /// Self type. I.e., `Vec<T,N>` if a vector, `Quat<T>` if a quaternion.
    typedef VType self_t;
    /// Element type.
    typedef T elem_t;
    /// Vector dimension.
    static constexpr index_t DIM = N;
    
    static const self_t ones;
    static const self_t zeros;
    static const self_t unit_x;

    /*******************************
     * Structors                   *
     *******************************/

protected:
    
    constexpr VecCommon():VecBase<T,N>() {}

    constexpr VecCommon(T a):VecBase<T,N>(a) {}

    constexpr VecCommon(const T a[N]):VecBase<T,N>(a) {}
    
    constexpr VecCommon(std::initializer_list<T> l):VecBase<T,N>(l.begin()) {}

    /*******************************
     * Operators                   *
     *******************************/

public:

    /**
     * @brief Vector element access.
     * 
     * @param idx Index of element to retrieve.
     * @return A read-only reference to the element at index `idx`.
     */
    inline const T& operator[](index_t idx) const {
        return this->get(idx);
    }
    
    /**
     * @brief Vector element access.
     * 
     * @param idx Index of element to retrieve.
     * @return A reference to the element at index `idx`.
     */
    inline T& operator[](index_t idx) {
        return this->get(idx);
    }

    /**
     * @brief Element typecast.
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

    /**
     * @brief Equality test.
     * @return `true` if all corresponding elements of `this` and `vv` are equal,
     *`false` otherwise.
     */
    inline bool operator==(const self_t& vv) const {
        for (index_t i = 0; i < N; i++) {
            if (vv.get(i) != this->get(i)) return false;
        }
        return true;
    }

#if __cplusplus < 202002L
    /**
     * @brief Inequality test.
     * @return `true` if any corresponding elements of `this` and `vv` are unequal,
     *`false` otherwise.
     */
    inline bool operator!=(const self_t& vv) const {
        return not ((*this) == vv);
    }
#endif

    /// Element-wise addition.
    inline self_t operator+(const self_t &v) const {
        return this->add(v);
    }

    /// Element-wise addition and assignment.
    self_t& operator+=(const self_t& vv) {
        for (index_t i = 0; i < N; i++) {
            this->get(i) += vv[i];
        }
        return *(static_cast<self_t*>(this));
    }
    
    /// Element-wise subtraction.
    inline self_t operator-(const self_t& v) const {
        return this->sub(v);
    }
    
    /// Subtraction and assignment.
    self_t& operator-=(const self_t& vv) {
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
    self_t& operator*=(const self_t& vv) {
        for (index_t i = 0; i < N; i++) {
            this->get(i) *= vv[i];
        }
        return *(static_cast<self_t*>(this));
    }

    /**
     * @brief Negation. 
     * @return A copy of this vector with all elements negated (i.e. a 
     * vector pointing in the opposite direction).
     */
    inline self_t operator-() const {
        self_t r;
        for (index_t i = 0; i < N; i++) {
            r[i] = -this->get(i);
        }
        return r;
    }

    /*******************************
     * Methods                     *
     *******************************/
    
    /**
     * @brief Vector addition.
     * @param v Another vector.
     * @return A new vector `x` such that `x[i] = this[i] + v[i]`.
     */
    inline self_t add(const self_t& v) const {
        self_t r;
        for (index_t i = 0; i < N; i++) {
            r[i] = this->get(i) + v.get(i);
        }
        return r;
    }
    
    /**
     * @brief Vector subtraction.
     * @param v Another vector.
     * @return A new vector `x` such that `x[i] = this[i] - v[i]`.
     */
    inline self_t sub(const self_t& v) const {
        self_t r;
        for (index_t i = 0; i < N; i++) {
            r[i] = this->get(i) - v[i];
        }
        return r;
    }
    
    /**
     * @brief Element-wise multiplication.
     * @param v Another vector.
     * @return A new vector `x` such that `x[i] = this[i] * v[i]`.
     */
    inline self_t scale(const self_t& v) const {
        self_t r;
        for (index_t i = 0; i < N; i++) {
            r[i] = this->get(i) * v[i];
        }
        return r;
    }
    
    /**
     * @brief Scalar multiple.
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
     * @brief Vector normalization.
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
     * @brief Dot product.
     * @param v Another vector.
     * @return The dot product of `this` with `v`.
     */
    inline T dot(const self_t& v) const {
        T sum = (*this)[0] * v[0];
        for (index_t i = 1; i < N; i++) {
            sum = std::fma(this->get(i), v[i], sum);
        }
        return sum;
    }
    
    /**
     * @brief Euclidean norm (magnitude).
     * @return The Euclidean magnitude (geometric length) of this vector.
     */
    inline T mag() const {
        return std::sqrt(this->mag2());
    }
    
    /**
     * @brief Squared magnitude.
     * @return The square of the magnitude of this vector.
     */
    inline T mag2() const {
        T x0  = this->get(0);
        T sum = x0 * x0;
        for (index_t i = 1; i < N; i++) {
            T xi = this->get(i);
            sum  = std::fma(xi, xi, sum);
        }
        return sum;
    }

    /**
     * @brief Distance between points.
     * @param pt Another point.
     * @return The distance between `this` and `pt`.
     */
    T dist(const self_t& pt) const {
        return (this->sub(pt)).mag();
    }

    /**
     * @brief Distance squared to a point.
     * @param pt Another point.
     * @return The square of the distance between `this` and `pt`.
     */
    T dist2(const self_t& pt) const {
        return (this->sub(pt)).mag2();
    }

    /**
     * @brief Reflection about a normal.
     * @param normal Axis of reflection.
     * @return A copy of this vector reflected across the given axis.
     */
    self_t reflect_about(self_t normal) const {
        normal = normal.unit();
        return this->neg() + normal * 2 * this->dot(normal);
    }

    /**
     * @brief Elastic collision.
     * 
     * Treat `this` as a velocity vector or incident ray; this function returns
     * the velocity reflected off of a surface with normal `normal`. 
     * Convenience for `-reflect_about(normal)`.
     * 
     * @param normal Normal of surface to "bounce" on.
     * @return The "bounced" direction vector.
     */
    self_t bounce_on(const self_t& normal) const {
        return -reflect_about(normal);
    }
    
    /**
     * @brief Orthogonal projection to an axis.
     * @param axis A direction vector.
     * @return A vector in direction `axis` with magnitude equal to the component
     * of `this` aligned with `axis`.
     */
    self_t project_on(const self_t& axis) const {
        T m2 = axis.mag2();
        return m2 != 0 ? (axis * this->dot(axis) / m2) : axis;
    }
    
    /**
     * @brief Return the component of `this` that projects to `axis`, as a fraction of axis's 
     * length.
     * @param axis An arbitrary basis vector.
     */
    T fraction_on(const self_t& axis) const {
        return this->dot(axis) / axis.mag2();
    }
    
    /**
     * @brief Compute a vector with the direction of `this` and a new magnitude `mag`.
     * 
     * If `this` is the zero vector, it will remain unchanged.
     */
    self_t with_length(T mag) const {
        T m2 = mag2();
        if (m2 == 0) return {};
        self_t r;
        for (index_t i = 0; i < N; ++i) {
            r[i] = mag * this->get(i) / std::sqrt(m2);
        }
        return r;
    }

    /**
     * @brief Linear interpolation. 
     * 
     * A mix parameter of 0 evaluates to `this`, while 1 is `v`.
     * @param v Another vector.
     * @param mix A mixing factor between 0 and 1.
     * @return A linear mixing of `this` with `v`.
     */
    self_t mix(const self_t& v, T mix) const {
        T xim = 1 - mix;
        Vec<T,N> r;
        for (index_t i = 0; i < N; i++) {
            r[i] = xim*this->get(i) + mix*v[i];
        }
        return r;
    }

    /**
     * @brief Angle between vectors.
     * @param v Another vector.
     * @return Angle in radians between `this` and `v`, between 0 and `pi`.
     */
    T angle_to(const self_t& v) const {
        T c = v.unit().dot(unit());
        c = std::min((T)1, std::max(c, (T)-1));
        return std::acos(c);
    }
    
    ///////////////////////
    // Clamping/Rounding //
    ///////////////////////
    
    /**
     * @brief Element-wise absolute value.
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
     * @brief Element-wise floor function.
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
     * @brief Element-wise ceiling function.
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
     * @brief Element-wise minimum of two `Vec`s.
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
     * @brief Element-wise maximum of two `Vec`s.
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
     * @brief Minimum element.
     * @return The value of the component with the lowest value.
     */
    T min() const {
        T x = (*this)[0];
        for (index_t i = 1; i < N; ++i) {
            x = std::min((*this)[i], x);
        }
        return x;
    }
    
    /**
     * @brief Maximum element.
     * @return The value of the component with the highest value.
     */
    T max() const {
        T x = (*this)[0];
        for (index_t i = 1; i < N; ++i) {
            x = std::max((*this)[i], x);
        }
        return x;
    }
    
    /**
     * @brief Element-wise clamp.
     * @param lo Element-wise lower extremes.
     * @param hi Element-wise upper extremes.
     * @return A new vector such that each element `x[i]` is clamped
     * between `lo[i]` and `hi[i]`.
     */
     self_t clamp(const self_t &lo, const self_t &hi) const {
        Vec<T,N> r;
        for (index_t i = 0; i < N; i++) {
            r[i] = std::min(std::max(this->get(i), lo[i]), hi[i]);
        }
        return r;
    }
    
    /**
     * @brief Round each element to the nearest integer.
     */
    self_t round() const {
        Vec<T,N> r;
        for (index_t i = 0; i < N; i++) {
            r[i] = std::round((*this)[i]);
        }
        return r;
    }
    
    /// Return the index of the coordinate with the largest absolute value.
    index_t argmax() const {
        T biggest = std::abs((*this)[0]);
        index_t k = 0;
        for (index_t i = 1; i < N; ++i) {
            T x = std::abs((*this)[i]);
            if (x > biggest) {
                biggest = x;
                k = i;
            }
        }
        return k;
    }
    
    /// Return the index of the coordinate with the smallest absolute value.
    index_t argmin() const {
        T smallest = std::abs((*this)[0]);
        index_t k = 0;
        for (index_t i = 1; i < N; ++i) {
            T x = std::abs((*this)[i]);
            if (x < smallest) {
                smallest = x;
                k = i;
            }
        }
        return k;
    }
    
    /**
     * @brief Resized copy of a vector.
     * @tparam M Dimension of new vector.
     * @return A new vector with size `M`. If `M` is larger than `N`, the
     * new elements will be set to zero. If `M` is 1, then the return type is `T`.
     */
    template <index_t M>
    inline typename PointType<T,M>::point_t resized() const {
        typename PointType<T,M>::point_t out;
        const T *p = this->begin();
        T* o = PointType<T,M>::iterator(out);
        std::copy(p, p + std::min(M,N), o);
        return out;
    }

    /**
     * @brief Return `true` if all elements are zero.
     */
    bool is_zero() const {
        for (index_t i = 0; i < N; i++) {
            if (this->get(i) != (T)0) return false;
        }
        return true;
    }

    /**
     * @brief The number of elements in this vector. Always equal to `N`. 
     */
    inline index_t size() const {
        return N;
    }

}; /* end VecCommon class */


template <typename T, index_t N, typename VType>
    const VType VecCommon< T,N,VType >::ones  = VType((T)1);
template <typename T, index_t N, typename VType>
    const VType VecCommon< T,N,VType >::zeros = VType((T)0);
template <typename T, index_t N, typename VType>
    const VType VecCommon< T,N,VType >::unit_x = []() {
        VType v((T)0);
        v[0] = (T)1;
        return v;
    }();

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
struct _RGBChannelConversion<
            T, 
            typename std::enable_if<std::is_integral<T>::value,void>::type
        >
{
    static const int scale = 1;
};

template <typename T, typename Enable=void>
struct IsVector {
    static const bool value = false;
};

template <typename T>
struct IsVector<T, typename std::enable_if<
                        std::is_base_of<
                            detail::VecCommon<typename T::elem_t, T::DIM, T>, 
                            T>::value,
                    void>::type
                 >
{
    const static bool value = true;
};


} /* end namespace detail */
} /* end namepsace geom */
