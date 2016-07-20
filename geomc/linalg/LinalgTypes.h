/*
 * LinalgTypes.h
 *
 *  This file shall not include any other files
 *
 *  Created on: Nov 11, 2010
 *      Author: tbabb
 */

#ifndef LINALGTYPES_H_
#define LINALGTYPES_H_

#include <geomc/geomc_defs.h>
#include <geomc/Storage.h>

/** @defgroup linalg Linalg
 *  @brief Linear algebra functions and classes.
 */

#include <boost/iterator/counting_iterator.hpp>

//vector convenience macros

/** @ingroup linalg 
 *  @{
 */

//4d
#define ONE_VEC4d   (Vec<double,4>::ones) 
#define ZERO_VEC4d  (Vec<double,4>::zeros)
#define X_AXIS4d    (Vec<double,4>::X_AXIS)
#define Y_AXIS4d    (Vec<double,4>::Y_AXIS)
#define Z_AXIS4d    (Vec<double,4>::Z_AXIS)
#define W_AXIS4d    (Vec<double,4>::W_AXIS)

#define ONE_VEC4f   (Vec<float,4>::ones)
#define ZERO_VEC4f  (Vec<float,4>::zeros)
#define X_AXIS4f    (Vec<float,4>::X_AXIS)
#define Y_AXIS4f    (Vec<float,4>::Y_AXIS)
#define Z_AXIS4f    (Vec<float,4>::Z_AXIS)
#define W_AXIS4f    (Vec<float,4>::W_AXIS)

//3d
#define ONE_VEC3d   (Vec<double,3>::ones)
#define ZERO_VEC3d  (Vec<double,3>::zeros)
#define X_AXIS3d    (Vec<double,3>::X_AXIS)
#define Y_AXIS3d    (Vec<double,3>::Y_AXIS)
#define Z_AXIS3d    (Vec<double,3>::Z_AXIS)

#define ONE_VEC3f   (Vec<float,3>::ones)
#define ZERO_VEC3f  (Vec<float,3>::zeros)
#define X_AXIS3f    (Vec<float,3>::X_AXIS)
#define Y_AXIS3f    (Vec<float,3>::Y_AXIS)
#define Z_AXIS3f    (Vec<float,3>::Z_AXIS)

//2d
#define ONE_VEC2d   (Vec<double,2>::ones)
#define ZERO_VEC2d  (Vec<double,2>::zeros)
#define X_AXIS2d    (Vec<double,2>::X_AXIS)
#define Y_AXIS2d    (Vec<double,2>::Y_AXIS)

#define ONE_VEC2f   (Vec<float,2>::ones)
#define ZERO_VEC2f  (Vec<float,2>::zeros)
#define X_AXIS2f    (Vec<float,2>::X_AXIS)
#define Y_AXIS2f    (Vec<float,2>::Y_AXIS)

// @}  // ingroup linalg

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
}

   typedef const void* storage_id_t;

    // data fwd decls
    template <typename T, index_t N>            class Vec;
    template <typename T>                       class Quat;
    template <typename T, index_t N>            class Ray;
    template <typename T, index_t N>            class AffineTransform;
    template <typename T, index_t M, index_t N> class PLUDecomposition;
    
    // mtx fwd decls
namespace detail {
    template <typename T, index_t M, index_t N, typename Derived> class MatrixBase;
}
    template <typename T, index_t M, index_t N, StoragePolicy P=STORAGE_UNIQUE>
                                                class SimpleMatrix;
    template <typename T, index_t M, index_t N> class DiagMatrix;
    template <typename Ma, typename Mb>         class AugmentedMatrix;
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
    typedef SimpleMatrix<double,4,4>  SimpleMatrix4d;
    typedef SimpleMatrix<double,3,3>  SimpleMatrix3d;
    typedef SimpleMatrix<double,2,2>  SimpleMatrix2d;
    typedef SimpleMatrix<double,0,0>  SimpleMatrixNd;
    
    typedef SimpleMatrix<float,4,4>   SimpleMatrix4f;
    typedef SimpleMatrix<float,3,3>   SimpleMatrix3f;
    typedef SimpleMatrix<float,2,2>   SimpleMatrix2f;
    typedef SimpleMatrix<float,0,0>   SimpleMatrixNf;
    
    // AffineTransform
    typedef AffineTransform<double,3> AffineTransform3d;
    typedef AffineTransform<double,2> AffineTransform2d;
    
    typedef AffineTransform<float,3>  AffineTransform3f;
    typedef AffineTransform<float,2>  AffineTransform2f;

#if __cplusplus >= 201103L

    template <typename T, index_t M, index_t N> 
    using WrapperMatrix = SimpleMatrix<T,M,N,STORAGE_USER_OWNED>;

#endif
    
    // Point type information
    // these are used to switch code between behaviors for
    // vectors and single-element types (float, int, ...)
    // this can be user-specialized, if need be, for interesting
    // types of T.
    
    // TODO: move me to the eventual matrix traits header.

    template <typename T, index_t N>
    struct PointType {
        typedef Vec<T,N> point_t;
        
        static inline T* iterator(point_t &p) {
            return p.begin();
        }
        
        static inline const T* iterator(const point_t &p) {
            return p.begin();
        }
        
        static inline point_t from_ptr(T* p) {
            return point_t(p);
        }
        
        static inline point_t from_larger_vector(const Vec<T,N+1> &v) {
            return v.template resized<N>();
        }
        
        static inline T mag2(const point_t &p) {
            return p.mag2();
        }
        
        static inline T mag(const point_t &p) {
            return p.mag();
        }
    };

    template <typename T>
    struct PointType<T,1> {
        typedef T point_t;
        
        static inline T* iterator(point_t &p) {
            return &p;
        }
        
        static inline const T* iterator(const point_t &p) {
            return &p;
        }
        
        static inline point_t from_ptr(T* p) {
            return *p;
        }
        
        static inline point_t from_larger_vector(const Vec<T,2> &v) {
            return v[0];
        }
        
        static inline T mag2(const point_t &p) {
            return p * p;
        }
        
        static inline T mag(const point_t &p) {
            return p;
        }
    };
    
} // namespace geom

#endif /* LINALGTYPES_H_ */
