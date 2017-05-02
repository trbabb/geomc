/*
 * Vec.h
 *
 *  Created on: Dec 25, 2010
 *      Author: tbabb
 */

#ifndef VEC_H_
#define VEC_H_

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_scalar.hpp>

#include <geomc/linalg/vecdetail/VecBase.h>  // the meat happens in here.
#include <geomc/linalg/vecdetail/Vec2.h>
#include <geomc/linalg/vecdetail/Vec3.h>
#include <geomc/linalg/vecdetail/Vec4.h>
#include <geomc/linalg/vecdetail/VecStd.h>

#include "mtxdetail/MatrixGlue.h"


namespace geom {
    
    // fwd decl
    namespace detail {
        template <typename M, typename RefType> class MtxColIterator;
    }

/*******************************
 * Vector Operators            *
 *******************************/

//TODO: convert to <auto>
//TODO: use boost::operator?

// why enable_if is necessary to prevent masking of
// more explicit matrix mult operators is beyond me.
// c++ is a dumb, awful language.

/**
 * Vector-scalar multiplication
 * 
 * @param v A vector
 * @param d Scalar value of type satisfying `boost::is_scalar`
 * @return A new vector `x` such that `x[i] = v[i] * d`
 * @related Vec
 */
#ifdef PARSING_DOXYGEN
    inline template <typename T, index_t N, typename U>
    Vec<T,N> operator* (const Vec<T,N> &v, U d) {}
#endif
template <typename V, typename U> 
inline typename boost::enable_if_c<
            boost::is_scalar<U>::value and detail::IsVector<V>::value,
        V>::type 
operator* (const V &v, U d) {
    V r;
    for (index_t i = 0; i < V::DIM; i++) {
        r[i] = v[i] * d;
    }
    return r;
}

/**
 * Vector-scalar multiplication
 * 
 * @param d Scalar value of type satisfying `boost::is_scalar`
 * @param v A vector
 * @return A new vector `x` such that `x[i] = d * v[i]`
 * @related Vec
 */
#ifdef PARSING_DOXYGEN
    inline template <typename T, index_t N, typename U>
    Vec<T,N> operator* (U d, const Vec<T,N> &v) {}
#endif
template <typename V, typename U> 
inline typename boost::enable_if_c<
        boost::is_scalar<U>::value and detail::IsVector<V>::value, 
        V>::type 
operator* (U d, const V &v) {
    V r;
    for (index_t i = 0; i < V::DIM; i++) {
        r[i] = d * v[i];
    }
    return r;
}

/**
 * Element-wise vector multiplication
 * 
 * @param a A vector
 * @param b A vector
 * @return A new vector `x` such that `x[i] = a[i] * b[i]`
 * @related Vec
 */
#ifdef PARSING_DOXYGEN
template <typename T, index_t N> 
const Vec<T,N> operator* (const Vec<T,N> &a, const Vec<T,N> &b) {}
#else
template <typename V>
inline
typename boost::enable_if_c<detail::IsVector<V>::value, V>::type 
operator*(const V &a, const V &b) {
    return a.scale(b);
}
#endif

/**
 * Vector division by a scalar
 * 
 * @param v A vector
 * @param d Scalar value
 * @return A new vector `x` such that `x[i] = v[i] / d`
 * @related Vec
 */
#ifdef PARSING_DOXYGEN
template <typename T, index_t N, typename U> 
inline Vec<T,N> operator/ (const Vec<T,N> &v, U d) {}
#else
template <typename V, typename U>
inline
typename boost::enable_if_c<detail::IsVector<V>::value, V>::type
operator/(const V &v, U d) {
    V r;
    for (index_t i = 0; i < V::DIM; i++) {
        r[i] = v[i] / d;
    }
    return r;
}
#endif

/**
 * Scalar division by a vector
 * 
 * @param d Scalar value
 * @param v A vector
 * @return A new vector `x` such that `x[i] = d / v[i]`
 * @related Vec
 */
#ifdef PARSING_DOXYGEN
template <typename T, index_t N, typename U> 
inline Vec<T,N> operator/ (const Vec<T,N> &v, U d) {}
#else
template <typename V, typename U>
inline
typename boost::enable_if_c<detail::IsVector<V>::value, V>::type
operator/(U d, const V &v) {
    V r;
    for (index_t i = 0; i < V::DIM; i++) {
        r[i] = d / v[i];
    }
    return r;
}
#endif

/**
 * Element-wise vector division
 * 
 * @param a A vector
 * @param b A vector
 * @return A new vector `x` such that `x[i] = a[i] / b[i]`
 * @related Vec
 */
#ifdef PARSING_DOXYGEN
template <typename T, index_t N> 
const Vec<T,N> operator/ (const Vec<T,N> &a, const Vec<T,N> &b) {}
#else
template <typename V>
typename boost::enable_if_c<detail::IsVector<V>::value,V>::type
operator/(const V &a, const V &b) {
    V r;
    for (index_t i = 0; i < V::DIM; i++) {
        r[i] = a[i] / b[i];
    }
    return r;
}
#endif

#ifdef GEOMC_LINALG_USE_STREAMS
template <typename T, index_t N> 
std::ostream &operator<< (std::ostream &stream, const Vec<T,N> &v) {
    stream << "(";
    for (index_t i = 0; i < N-1; i++) {
        stream << v[i] << ", ";
    }
    stream << v[N-1] << ")";
    return stream;
}

template <typename T>
std::ostream &operator<<(std::ostream &stream, const Quat<T> &q) {
    stream << "quat(";
    stream << q.x << ", ";
    stream << q.y << ", ";
    stream << q.z << " | ";
    stream << q.w << ")";
    return stream;
}
#endif


/* Declare the Vec type that we'll actually use, as a subclass of the hidden base template class.
 * We do this so that all sizes of vector can share common code, since a template that is merely specialized
 * does not implicitly inherit functionality from the base template. With this method, Vec<T,2> can be a specialization
 * of Vec<T,N>, but still keep the functionality of VecBase<T,N>.
 * 
 * All vectors are sizeof(T)*N and can be packed.
 */

/**
 * @ingroup linalg
 * @brief A tuple of `N` elements of type `T`.
 * 
 * Vectors are lightweight and generally perform as well as a bare array of their
 * element type. 
 * 
 * Geomc makes no type distinction between vectors, points, or normals; the 
 * distinction is to be made by the programmer based on usage.
 * 
 * Declaring a 3-dimensional vector of doubles:
 * 
 *     Vec<double,3> v;
 * 
 * Basic arithmetic:
 *     
 *     v3 = v1 + v2;   // addition
 *     v3 = v1 - v2;   // subtraction
 *     v3 = -v1;       // negation
 *     v3 = 2.71 * v1; // scalar mult
 *     v3 = 1.61 / v1; // scalar div
 *     v3 = v1 / 1.41; // scalar div
 * 
 * Element-wise multiplication / division:
 * 
 *     v3 = v1 * v2;
 *     v3 = v2 / v1;
 * 
 * Element access:
 * 
 *     double z = v[2];
 *     z = v.get(2);
 * 
 * Access to / iteration over internal array:
 * 
 *     for (double *p = v.begin(); p != v.end(); p++) {
 *         double a = f(*p, ...);
 *     }
 * 
 * Cross product (3D only):
 * 
 *     v3 = v1 ^ v2;
 * 
 * Compatibility with std:
 * 
 *     // element-wise operations
 *     std::min(v1, v2);
 *     std::max(v1, v2);
 *     std::abs(v1);
 *     std::floor(v1);
 *     std::ceil(v1);
 * 
 * Resizing:
 *     
 *     Vec<double,3> v3d;
 *     Vec<double,2> v2d = v3d.resized<2>(); // truncate the last coordinate
 *     Vec<double,4> v4d = v3d.resized<4>(); // last coordinate is zero
 * 
 */
template <typename T, index_t N> class Vec : public detail::VecCommon< T,N,Vec<T,N> > {
public:
    
    /**
     * Construct a new vector with all elements set to zero.
     */
    Vec():detail::VecCommon< T,N,Vec<T,N> >() {}

    /**
     * Construct a new vector with all elements set to the value of `a`.
     * 
     * @param a Scalar value
     */
    Vec(T a):detail::VecCommon< T,N,Vec<T,N> >(a) {}
    
    /**
     * Construct a new vector with elements copied from `a`.
     * @param a An array of length `N`. 
     */
    Vec(const T a[N]):detail::VecCommon< T,N,Vec<T,N> >(a) {}
    
    /**
     * Construct a new vector with the elements from `v`, with `a` as the last 
     * element.
     * @param v A vector of dimension `N - 1`
     * @param a The value of the last element
     */
    Vec(const Vec<T,N-1> &v, T a) {
        std::copy(v.begin(), v.end(), detail::VecBase<T,N>::begin());
        this->get(N-1) = a;
    }
    
    /**
     * Construct a vector from a column of a matrix.
     * @param mtx_col A matrix column iterator (obtained via `mtx.col(i)`).
     */
    template <typename Mx, typename Ref>
    Vec(detail::MtxColIterator<Mx,Ref> mtx_col) {
        T *p = this->begin();
        const index_t n = std::min(N, mtx_col->mtx->cols());
        for (index_t i = 0; i < n; i++, p++, mtx_col++) {
            *p = *mtx_col;
        }
    }
    
    /**
     * Construct a vector from a brace-initialization list. (c++11)
     *
     * Example: `Vec<int,3> v = {2, 5, 8};`
     * @param items A brace-initializer list.
     */
#if __cplusplus >= 201103L or PARSING_DOXYGEN
    Vec(const std::initializer_list<T>& items):detail::VecCommon< T,N,Vec<T,N> >(items.begin()) {
#if __cplusplus >= 201402L
        // items.size() not constexpr in c++11  D:<
        static_assert(items.size() == N);
#endif
    }
#endif
    
}; /* class Vec */

} // namespace geom


#endif /* VEC_H_ */
