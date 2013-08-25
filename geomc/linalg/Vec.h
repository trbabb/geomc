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


namespace geom {
    
    
/** @ingroup linalg
 *  @{
 */

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
 */
#ifdef PARSING_DOXYGEN
    inline template <typename T, index_t N, typename U>
    Vec<T,N> operator* (const Vec<T,N> &v, U d) {}
#endif
template <typename T, index_t N, typename U> 
inline typename boost::enable_if<boost::is_scalar<U>,Vec<T,N> >::type operator* (const Vec<T,N> &v, U d) {
    Vec<T,N> r;
    for (index_t i = 0; i < N; i++){
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
 */
#ifdef PARSING_DOXYGEN
    inline template <typename T, index_t N, typename U>
    Vec<T,N> operator* (U d, const Vec<T,N> &v) {}
#endif
template <typename T, index_t N, typename U> 
inline typename boost::enable_if<boost::is_scalar<U>,Vec<T,N> >::type operator* (U d, const Vec<T,N> &v) {
    Vec<T,N> r;
    for (index_t i = 0; i < N; i++){
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
 */
template <typename T, index_t N> 
const Vec<T,N> operator* (const Vec<T,N> &a, const Vec<T,N> &b) {
    return a.scale(b);
}

/**
 * Vector division by a scalar
 * 
 * @param v A vector
 * @param d Scalar value
 * @return A new vector `x` such that `x[i] = v[i] / d`
 */
template <typename T, index_t N, typename U> 
inline Vec<T,N> operator/ (const Vec<T,N> &v, U d) {
    Vec<T,N> r;
    for (index_t i = 0; i < N; i++){
        r[i] = v[i] / d;
    }
    return r;
}

/**
 * Scalar division by a vector
 * 
 * @param d Scalar value
 * @param v A vector
 * @return A new vector `x` such that `x[i] = d / v[i]`
 */
template <typename T, index_t N, typename U> 
inline Vec<T,N> operator/ (U d, const Vec<T,N> &v) {
    Vec<T,N> r;
    for (index_t i = 0; i < N; i++){
        r[i] = d / v[i];
    }
    return r;
}

/**
 * Element-wise vector division
 * 
 * @param a A vector
 * @param b A vector
 * @return A new vector `x` such that `x[i] = a[i] / b[i]`
 */
template <typename T, index_t N> 
const Vec<T,N> operator/ (const Vec<T,N> &a, const Vec<T,N> &b) {
    Vec<T,N> r;
    for (index_t i = 0; i < N; i++){
        r[i] = a[i] / b[i];
    }
    return r;
}

#ifdef GEOMC_LINALG_USE_STREAMS
template <typename T, index_t N> 
std::ostream &operator<< (std::ostream &stream, const Vec<T,N> &v) {
    stream << "(";
    for (index_t i = 0; i < N-1; i++){
        stream << v[i] << ", ";
    }
    stream << v[N-1] << ")";
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
 */
template <typename T, index_t N> class Vec : public detail::VecCommon<T,N> {
public:
    
    /**
     * Construct a new vector with all elements set to zero.
     */
    Vec():detail::VecCommon<T,N>(){}

    /**
     * Construct a new vector with all elements set to the value of `a`.
     * 
     * @param a Scalar value
     */
    Vec(T a):detail::VecCommon<T,N>(a){}
    
    /**
     * Construct a new vector with elements copied from `a`.
     * @param a An array of length `N`. 
     */
    Vec(T a[N]):detail::VecCommon<T,N>(a){}
    
    /**
     * Construct a new vector with the elements from `v`, with `a` as the last 
     * element.
     * @param v A vector of dimension `N - 1`
     * @param a The value of the last element
     */
    Vec(Vec<T,N-1> &v, T a){
        std::copy(v.begin(), v.end(), detail::VecBase<T,N>::begin());
        this->get(N-1) = a;
    }
    
}; /* class Vec */

}; // namespace geom


namespace std {

// these make vectors more interchangeable with bare types:
/**
 * Elemnt-wise maximum.
 */
template <typename T, index_t N>
inline geom::Vec<T,N> max(const geom::Vec<T,N> &a, const geom::Vec<T,N> &b) {
    return a.max(b);
}

/**
 * Element-wise minimum.
 */
template <typename T, index_t N>
inline geom::Vec<T,N> min(const geom::Vec<T,N> &a, const geom::Vec<T,N> &b) {
    return a.min(b);
}

/**
 * Element-wise absolute value.
 */
template <typename T, index_t N>
inline geom::Vec<T,N> abs(const geom::Vec<T,N> &v) {
    return v.abs();
}

/**
 * Element-wise floor.
 */
template <typename T, index_t N>
inline geom::Vec<T,N> floor(const geom::Vec<T,N> &v) {
    return v.floor();
}

/**
 * Element-wise ceiling.
 */
template <typename T, index_t N>
inline geom::Vec<T,N> ceil(const geom::Vec<T,N> &v) {
    return v.ceil();
}

}; // namespace std

/// @} //ingroup linalg


#endif /* VEC_H_ */
