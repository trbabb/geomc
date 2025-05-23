#pragma once
/*
 * LinalgTypes.h
 *
 *  Created on: Nov 11, 2010
 *      Author: tbabb
 */

#include <geomc/geomc_defs.h>
#include <geomc/Storage.h>

/** @defgroup linalg Linalg
 *  @brief Linear algebra functions and classes.
 */

// type definitions

namespace geom {

namespace detail {

enum VecOrientation {
    ORIENT_VEC_ROW,
    ORIENT_VEC_COL,
    ORIENT_VEC_UNKNOWN
};

enum MatrixResultAgreement {
    MTX_RESULT_MATCH,
    MTX_RESULT_MISMATCH,
    MTX_RESULT_UNKNOWN
};

} // namespace detail

enum MatrixLayout {
    ROW_MAJOR,
    COL_MAJOR
};

// represents a memory range:
struct storage_token_t {
    const void* lo;
    const void* hi;
};

// data fwd decls
template <typename T, index_t N>            class Vec;
template <typename T>                       class Quat;
template <typename T, index_t N>            class Ray;
template <typename T, index_t N>            class Rotation;
template <typename T, index_t N>            class Isometry;
template <typename T, index_t N>            class Similarity;
template <typename T, index_t N>            class AffineTransform;
template <typename T, index_t M, index_t N> class PLUDecomposition;

// mtx fwd decls
namespace detail {
    template <typename T, index_t M, index_t N, typename Derived> class MatrixBase;
}
template <
    typename T, 
    index_t M, 
    index_t N,
    MatrixLayout Lyt=ROW_MAJOR, 
    StoragePolicy P=STORAGE_UNIQUE
    >
class SimpleMatrix;
template <typename T, index_t M, index_t N> class DiagMatrix;
template <index_t N>                        class PermutationMatrix;
template <typename T>                       class MatrixHandle;

// Vec
typedef Vec<double,2>  Vec2d;
typedef Vec<float,2>   Vec2f;
typedef Vec<index_t,2> Vec2i;

typedef Vec<double,3>  Vec3d;
typedef Vec<float,3>   Vec3f;
typedef Vec<index_t,3> Vec3i;

typedef Vec<double,4>  Vec4d;
typedef Vec<float,4>   Vec4f;
typedef Vec<index_t,4> Vec4i;

typedef Quat<double>   Quatd;
typedef Quat<float>    Quatf;

// Ray
typedef Ray<double,4>  Ray4d;
typedef Ray<double,3>  Ray3d;
typedef Ray<double,2>  Ray2d;

typedef Ray<float,4>   Ray4f;
typedef Ray<float,3>   Ray3f;
typedef Ray<float,2>   Ray2f;

typedef Ray<index_t,4> Ray4i;
typedef Ray<index_t,3> Ray3i;
typedef Ray<index_t,2> Ray2i;

// Matrix
typedef SimpleMatrix<double,4,4> SimpleMatrix4d;
typedef SimpleMatrix<double,3,3> SimpleMatrix3d;
typedef SimpleMatrix<double,2,2> SimpleMatrix2d;
typedef SimpleMatrix<double,0,0> SimpleMatrixNd;

typedef SimpleMatrix<float,4,4>  SimpleMatrix4f;
typedef SimpleMatrix<float,3,3>  SimpleMatrix3f;
typedef SimpleMatrix<float,2,2>  SimpleMatrix2f;
typedef SimpleMatrix<float,0,0>  SimpleMatrixNf;

// AffineTransform
typedef AffineTransform<double,3> AffineTransform3d;
typedef AffineTransform<double,2> AffineTransform2d;

typedef AffineTransform<float,3>  AffineTransform3f;
typedef AffineTransform<float,2>  AffineTransform2f;

template <typename T, index_t M, index_t N, MatrixLayout Lyt=ROW_MAJOR> 
using WrapperMatrix = SimpleMatrix<T,M,N,Lyt,STORAGE_WRAPPED>;
    
/// Represents a means of transforming a point or vector.
template <typename Xf, typename T, index_t N>
concept Transform = requires (Xf xf, Vec<T,N> p) {
    // transforms can be applied to points
    { xf * p } -> std::convertible_to<Vec<T,N>>;
    // inverse transforms can be applied to points
    { p / xf } -> std::convertible_to<Vec<T,N>>;
    // transforms can be composed
    { xf * xf } -> std::convertible_to<Xf>;
    // transforms have inverses
    { xf.inverse() } -> std::convertible_to<Xf>;
    // transforms can operate on directions
    { xf.apply_direction(p) } -> std::convertible_to<Vec<T,N>>;
    // inverse transforms can operate on directions
    { xf.apply_inverse_direction(p) } -> std::convertible_to<Vec<T,N>>;
    // transforms can in-place compose
    { xf *= xf } -> std::convertible_to<Xf>;
    // transforms can in-place apply inverse
    { xf /= xf } -> std::convertible_to<Xf>;
};

/// A transformable is that which can be operated on by a Transform
/// and maintain its type.
template <typename S, typename Xf>
concept Transformable = 
    Transform<Xf, typename S::elem_t, S::N> and
    requires (S s, Xf xf) {
        // transformables can be transformed
        { xf * s } -> std::convertible_to<S>;
        // transformables can be inverse-transformed
        { s / xf } -> std::convertible_to<S>;
    };

// Point type information
// these are used to switch code between behaviors for
// vectors and single-element types (float, int, ...)
// this can be user-specialized, if need be, for interesting
// types of T.

// todo: should I live in Vec.h?

template <typename T, index_t N>
struct PointType {
    using point_t = Vec<T,N>;
    
    static inline T* iterator(point_t& p) {
        return p.begin();
    }
    
    static inline const T* iterator(const point_t& p) {
        return p.begin();
    }
    
    static inline point_t from_ptr(T* p) {
        return point_t(p);
    }
    
    static inline point_t from_larger_vector(const Vec<T,N+1>& v) {
        return v.template resized<N>();
    }
    
    static inline T mag2(const point_t& p) {
        return p.mag2();
    }
    
    static inline T mag(const point_t& p) {
        return p.mag();
    }
    
    static inline T dot(const point_t& p0, const point_t& p1) {
        return p0.dot(p1);
    }
    
    static inline point_t unit(const point_t& p) {
        return p.unit();
    }
};

template <typename T>
struct PointType<T,1> {
    using point_t = T;
    
    static inline T* iterator(point_t& p) {
        return &p;
    }
    
    static inline const T* iterator(const point_t& p) {
        return &p;
    }
    
    static inline point_t from_ptr(T* p) {
        return *p;
    }
    
    static inline point_t from_larger_vector(const Vec<T,2>& v) {
        return v[0];
    }
    
    static inline T mag2(const point_t& p) {
        return p * p;
    }
    
    static inline T mag(const point_t& p) {
        return std::abs(p);
    }
    
    static inline T dot(const point_t& p0, const point_t& p1) {
        return p0 * p1;
    }
    
    static inline T unit(const point_t& p) {
        return p >= 0 ? (T)1 : (T)-1;
    }
};

/// @brief Represents an object or operation that exists in a certain
/// dimension with a certain coordinate type.
/// @ingroup linalg
template <typename T, index_t _N>
struct Dimensional {
    /// The coordinate type of this object.
    using elem_t = T;
    /// The dimension of this object.
    static constexpr index_t N = _N;
    /**
        * @brief The type of a point in this object's space.
        *
        * An N-vector of T if N > 1, otherwise a T.
        */
    using point_t = PointType<T,N>::point_t;
    
private:
    // prevent this class from taking space
    int _dummy[0];
};

/**
    * @brief The type of a vector in N dimensions with elements of type T.
    *
    * VecType<T,N> is a Vec<T,N> in dimension > 1, otherwise a T.
    */
template <typename T, index_t N>
using VecType = typename PointType<T,N>::point_t;

/**
    * @brief Access the `n`th coordinate of a vector or scalar.
    */
template <typename T, index_t N>
inline constexpr T& coord(Vec<T,N>& v, index_t n) {
    return v[n];
}

/**
    * @brief Access the `n`th coordinate of a const vector or scalar.
    */
template <typename T, index_t N>
inline constexpr T coord(const Vec<T,N>& v, index_t n) {
    return v[n];
}

/**
    * @brief Access the `n`th coordinate of a 1-D vector (scalar).
    */
template <typename T>
inline constexpr T& coord(T& v, index_t n) {
    return v;
}

/**
    * @brief Access the `n`th coordinate of a const 1-D vector (scalar).
    */
template <typename T>
inline constexpr T coord(const T& v, index_t n) {
    return v;
}

/**
    * @brief Compute the squared magnitude of a vector or scalar.
    */
template <typename T, index_t N>
inline constexpr T mag2(const Vec<T,N>& v) {
    return v.mag2();
}

/**
    * @brief Compute the square of a scalar.
    */
template <typename T>
inline constexpr T mag2(T x) {
    return x * x;
}

/**
    * @brief Compute the magnitude of a vector or scalar.
    */
template <typename T, index_t N>
inline constexpr T mag(const Vec<T,N>& v) {
    return v.mag();
}

/**
    * @brief Compute the absolute value of a scalar.
    */
template <typename T>
inline constexpr T mag(T x) {
    return std::abs(x);
}

/**
    * @brief Concept for an object or operation that exists in a certain
    * dimension with a certain coordinate type.
    *
    * Inheriting from Dimensional is sufficient to satisfy this concept.
    */
template <typename Obj>
concept DimensionalObject = (Obj::N > 0) and requires {
    typename Obj::elem_t;
    { Obj::N } -> std::convertible_to<index_t>;
    typename Obj::point_t;
};

template <typename Obj, typename T, index_t N>
concept NDimensional = DimensionalObject<Obj> and
    std::same_as<typename Obj::elem_t, T> and
    Obj::N == N;

template <typename T>
struct VectorDimension {
    static constexpr index_t N = 1;
};

template <typename T, index_t M>
struct VectorDimension<Vec<T,M>> {
    static constexpr index_t N = M;
};

namespace detail {

template <typename T>
struct CoordType {
    using elem_t = T; // default to T for scalar types
};

template <DimensionalObject T>
struct CoordType<T> {
    // for dimensional objects, use the elem_t defined in the object
    using elem_t = typename T::elem_t;
};

} // namespace detail

template <typename T>
using Coordinate = typename detail::CoordType<T>::elem_t;

} // namespace geom
