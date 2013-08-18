/* Matrix.h
 *  
 * Created on: Oct 27, 2010
 *      Author: tbabb
 */

/** @ingroup linalg
 *  @defgroup matrix Matrix
 * 
 * @brief Matrix-related functions and classes
 * 
 * Include 
 * =======
 * `#include <geomc/linalg/Matrix.h>`
 * 
 * Design
 * ======
 * There are currently six distinct matrix template classes.
 * They are all interoperable, and provide iterators which
 * are functionally interchangeable with each other and
 * with pointers (row iterators, column iterators, region
 * iterators, and row-major matrix body iterators).
 * 
 * These iterators are compatible with std::copy(). Some
 * of them are not writeable (or may possibly throw an
 * error upon writing) in the case where a matrix element
 * does not have a corresponding memory location (DiagonalMatrix
 * is one such case, for example; off-diagonals are not stored).
 * 
 * In general, this scheme was designed to satisfy the following
 * requirements:
 * 
 * - All matrix types must interoperate.
 * - Using matrices with arithmetic operators should be straightforward,
 *   efficient, and feel like using a native type.
 * - Dimension mismatches should be caught at compile time wherever possible,
 *   to avoid unnecessary runtime checks.
 * - Element access shall be non-virtual and inline-able wherever possible.
 * - Dynamic memory allocations should be minimized or eliminated wherever possible.
 * - Copies to and from contiguous memory shall be fast where the
 *   internal matrix representation is contiguous.
 * - Handles to matrices of arbitrary type are possible.
 * 
 * Use
 * ===
 * 
 * Matrix dimensions
 * -----------------
 * 
 * Most matrices have templated size:
 * 
 *     SimpleMatrix<double, 3, 4> mat3x4;
 * 
 * In the example above, we construct a matrix with 3 rows and 4 columns. This 
 * matrix has **static dimensions**, in that its size is chosen (and fixed) at 
 * compile-time. Matrices with dimensions chosen at runtime may be written as:
 * 
 *     SimpleMatrix<double, 0, 0> matNxN(nrows, ncols);
 * 
 * There are two important functional differences between static matrices and 
 * dynamic matrices:
 * 
 *   - Static matrix operators may perform dimension "agreement" checks at 
 *     compile-time, catching errors early, and circumventing the cost of runtime 
 *     dimension checks.
 *   - Static matrices use "copy" semantics, while dynamic matrices use 
 *     "reference" semantics, and point to common storage when copy-constructed.
 * 
 * The latter point is important to understand, as writing:
 * 
 *     SimpleMatrix<double, 0, 0> m1 = some_matrix_function(...);
 *     SimpleMatrix<double, 0, 0> m2 = m1;
 *     m2[1][2] = 123;
 * 
 * means that element `(1, 2)` is altered in *both* `m1` and `m2`. This is **not** 
 * the case with static matrices.
 * 
 * To check if two matrices share common storage:
 *     
 *     if (mtx_aliases_storage(m1, m2)) { ... }
 * 
 * To copy the elements from one matrix to another, without referencing, and regardless
 * of matrix type or static/dynamic setting:
 * 
 *     mtxcopy(&dst_mtx, src_mtx);
 * 
 * Because of difference in storage strategy between static and dynamic matrices, copy 
 * and pass-by-value operations can be much heavier for large static matrices, 
 * because the entire matrix must be copied. Therefore it may be advisable to 
 * choose dynamic matrices for any data much larger than a few elements along 
 * each axis.
 * 
 * Dimension checking
 * ------------------
 * 
 * Matrix operations often place restrictions on the dimensions of their matrix 
 * operands. For example, matrix addition requires that both operands have the 
 * same dimensions (in other words, a 3 x 4 matrix can only be added to another 
 * 3 x 4 matrix). 
 * 
 * When matrices with static dimensions are involved, the library can often prove
 * at compile time that an argument's dimensions are correct or incorrect. If proven
 * correct, run-time dimension checks can be skipped, resulting in slightly faster
 * code. If incorrect, the compiler will error, catching program correctness problems
 * early. Compile-time dimension mismatches / requirement failures will generally
 * manifest as `"template argument deduction/substitution failed"` errors.
 * 
 * This code will prove its dimension-correctness at compile time:
 * 
 *     SimpleMatrix<double, 4, 4> mat4x4;
 *     DiagMatrix<double, 4, 4> dmat4x4;
 *     mat4x4 + dmat4x4  // statically proven; no runtime dimension checking
 * 
 * These will defer to a runtime check, since the dimensions of some arguments
 * cannot be deduced from their type:
 * 
 *     SimpleMatrix<double, 4, 4> mat4x4;
 *     SimpleMatrix<double, 0, 0> mat_a_NxN(4,4); // 4x4 matrix
 *     SimpleMatrix<double, 0, 0> mat_b_NxN(3,2); // 3x2 matrix
 *     
 *     mat4x4 + mat_a_NxN; // runtime dimension check; will succeed.
 *     mat4x4 + mat_b_NxN; // runtime dimension check; will throw an exception.
 * 
 * 
 * Runtime dimension checks will throw either a `DimensionMismatchException` or
 * `NonsquareMatrixException` on failure. Runtime checks can be disabled
 * completely (at the hazard of introducing memory access violations and other bugs)
 * by un-defining `GEOMC_MTX_CHECK_DIMS` in `geomc_defs.h`. A runtime dimension check
 * will occurr if any checked dimension is dynamic.
 * 
 * Operators
 * ---------
 * 
 * Matrices support most basic arithmetic operators. We'll demonstrate with these 
 * example objects:
 * 
 *     SimpleMatrix<double,3,3> m1;
 *     SimpleMatrix<double,3,4> m2;
 *     Vec<double,3> v;
 * 
 * Inter-matrix mult:
 * 
 *     m1 * m2
 * 
 * Matrix-scalar mult:
 * 
 *     m1 * 3
 *     1.618 * m1
 * 
 * Matrix-vector mult:
 * 
 *     m1 * v
 *     v * m2
 * 
 * Matrix addition / subtraction:
 * 
 *     m1 + m1
 *     m1 - m1
 * 
 * Equality test:
 * 
 *     m1 == m1
 *     m1 != m2
 * 
 * Indexing:
 * 
 *     double x = m1[1][2];
 *     m[0][1] = 2.718;
 * 
 * Accessing elements
 * ------------------
 * 
 * Indexing:
 * 
 *     // equivalent:
 *     float f1 = m.get(2, 3);
 *     float f2 = m[2][3]; 
 * 
 * Assignment:
 *     
 *     // all equivalent:    
 *     m.set(2, 3, val);
 *     m[2][3] = val;
 *     m.get(2, 3) = val;
 * 
 * Matrix body iterators:
 *     
 *     typedef SimpleMatrix<double,3,3> mat3;
 *     mat3 m;
 *     // iterate over the matrix body in row-major order:
 *     for (mat3::iterator i = mat.begin(), i != m.end(); m++) {
 *         *i = ... ;
 *     }
 * 
 * Matrix region iterators:
 * 
 *     typedef SimpleMatrix<double,3,3> mat3;
 *     mat3 m;
 *     Rect2i r = Rect2i(Vec2i(1,1), Vec2i(3,2));
 *     // iterate over the elements in region `r` in row-major order:
 *     for (mat3::region_iterator i = m.region_begin(r); i != m.region_end(r); i++) {
 *         Vec2i c = i.point();
 *         *i = f(c, ...);
 *     }
 * 
 * .
 */

//TODO: make exposed templates for the dynamically-chosen return types.
//TODO: add parameters to all the docs
//TODO: document:
//      x matrix static/dynamic dims
//      x matrix allowed operators
//      x document inter-operation / dim checking
//TODO: error reporting templates
//      x not fixing. c++ is too broken to support this (extremely basic) feature.
//TODO: rename Matrix (and SimpleMatrix?)
//TODO: rename MatrixRegion to vec2i, rect2i, etc.
//TODO: create a 'subregion' matrix.
//TODO: change 'scale' to 'mul'.
//TODO: get 'intersection/intersects/overlap/union/etc' terminology in Rect consistent.


//TODO: templatize memory layout choice?
//TODO: reduce bloat, particularly in matrix inv case.
//TODO: determinant
//TODO: matrix kernel (ND normal-finder)
//TODO: clean arbitrary Matrix handle construction
//TODO: verify correct matrix template function resolution.

//future: look for ways to keep object code bloat down.
//future: mechanism for col-major mtxs? (iterator based wrapper?) 
//        (could clutter code, because then other classes must be similarly templated).

//refactoring:

//TODO: iterators for nonzero entries.
//TODO: set() should return reference to (or value of) new element, so that `z = (mtx[x][y] = foo)` can
//      evaluate properly with proxy references.

#ifndef MATRIX_H_
#define MATRIX_H_

#include <boost/type_traits/is_scalar.hpp>

#include <geomc/linalg/LinalgTypes.h>

// matrix types
// each includes matrixdetail and matrixglue

#include <geomc/linalg/mtxtypes/SimpleMatrix.h>
#include <geomc/linalg/mtxtypes/AugmentedMatrix.h>
#include <geomc/linalg/mtxtypes/DiagMatrix.h>
#include <geomc/linalg/mtxtypes/SparseMatrix.h>
#include <geomc/linalg/mtxtypes/MatrixHandle.h>
#include <geomc/linalg/mtxtypes/PermutationMatrix.h>

#include <geomc/linalg/mtxdetail/MatrixFunctionImpl.h>
#include <geomc/linalg/mtxdetail/MatrixMult.h>
#include <geomc/linalg/mtxdetail/MatrixInv.h>
#include <geomc/linalg/mtxdetail/MatrixArithmetic.h>

#ifdef GEOMC_LINALG_USE_STREAMS
  #include <iostream>
  #include <iomanip>

#include "mtxdetail/MatrixGlue.h"
#endif


namespace geom {
    
/** @addtogroup matrix
 *  @{
 */
    
/************************************
 * Storage aliasing check           *
 ************************************/

/**
 * Do two matrices / vectors share storage?
 * 
 * In other words, might writing to one object change the contents of the other?
 * `Ma` and `Mb` must be matrix or vector types.
 * 
 * @param [in] a A matrix or vector object.
 * @param [in] b A matrix or vector object.
 * 
 * @return `true` if writing to `a` may alter `b` or vice versa; `false` otherwise.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> bool mtx_aliases_storage(const Ma &a, const Mb &b) {}
#endif
// do two matrices share any storage?
template <typename Ma, typename Mb>
typename boost::enable_if_c<
            detail::IsMatrix<Ma>::val and detail::IsMatrix<Mb>::val, 
            bool>::type 
mtx_aliases_storage(const Ma &a, const Mb &b) {
    // todo: currently assumes no /internal/ aliasing.
    typename Ma::storagebuffer_t buf_a = a.getStorageIDBuffer();
    typename Mb::storagebuffer_t buf_b = b.getStorageIDBuffer();
    
    a.getStorageIDs(buf_a.get()); // may be dynamically allocated
    b.getStorageIDs(buf_b.get()); // but usually static, and usually only one ptr
    for (index_t i = 0; i < a.getStorageIDCount(); i++) {
        for (index_t j = 0; j < b.getStorageIDCount(); j++) {
            if (buf_a.get()[i] == buf_b.get()[j]) return true;
        }
    }
    
    return false;
}

// needed for checks with vectors
template <typename Ma, typename Mb>
typename boost::enable_if_c<
            not (detail::IsMatrix<Ma>::val and detail::IsMatrix<Mb>::val), 
            bool>::type 
inline mtx_aliases_storage(const Ma &a, const Mb &b) {
    return false;
}

template <typename T, index_t N>
inline bool mtx_aliases_storage(const Vec<T,N> &a, const Vec<T,N> &b) {
    return a.begin() == b.begin();
}

/*********************************
 * Matrix copy                   *
 *********************************/

/** 
 * Copy the contents of `src` into `into`. This function may be optimized to 
 * perform better than running `std::copy()` on some matrix types' iterators. It
 * will also succeed in certain situations where `std::copy()` would fail (for
 * example when copying the contents of a diagonal matrix to another, any attempt
 * to write off the diagonal will throw an exception).
 * 
 * `Md` and `Ms` must be matrix or vector types whose dimensions match.
 * If the dimensions can be determined to mismatch at compile-time, the program
 * is considered invalid and the compilation will fail. If either object has dynamic
 * size, the check will be deferred to runtime, throwing a `DimensionMismatchException`
 * if the check fails.
 * 
 * @param [out] into A writeable matrix with dimensions matching `src`'s.
 * @param [in]  src  A matrix object.
 */
#ifdef PARSING_DOXYGEN
template <typename Md, typename Ms> void mtxcopy(Md *into, const Ms &src) {}
#endif
template <typename Md, typename Ms>
void mtxcopy(Md *into, const Ms &src,
                typename boost::enable_if_c<
                     detail::LinalgDimensionMatch<Md,Ms>::val and
                     (detail::_ImplVecOrient<Md,Ms>::orient != detail::ORIENT_VEC_UNKNOWN), 
                     int>::type dummy=0) {
#ifdef GEOMC_MTX_CHECK_DIMS
    typedef detail::_ImplMtxAdaptor<Md, detail::_ImplVecOrient<Md,Ms>::orient> D;
    typedef detail::_ImplMtxAdaptor<Ms, detail::_ImplVecOrient<Md,Ms>::orient> S;
    
    // runtime dimension match check
    if ((D::ROWDIM * S::ROWDIM == DYNAMIC_DIM and D::rows(*into) != S::rows(*into)) or
        (D::COLDIM * S::COLDIM == DYNAMIC_DIM and D::cols(*into) != S::cols(*into))) {
        throw DimensionMismatchException(D::rows(*into), D::cols(*into), S::rows(src), S::cols(src));
    }
#endif
    // todo: check mem aliasing?
    detail::_mtxcopy(into, src);
}


// unknown vector orientation case (dynamic matrix <-> vector) 
template <typename Md, typename Ms>
void mtxcopy(Md *into, const Ms &src,
                 typename boost::enable_if_c<
                    detail::_ImplVecOrient<Md,Ms>::orient == detail::ORIENT_VEC_UNKNOWN, 
                 int>::type dummy=0) {
#ifdef GEOMC_MTX_CHECK_DIMS
    typedef detail::_ImplMtxAdaptor<Md, detail::ORIENT_VEC_COL> Dc;
    typedef detail::_ImplMtxAdaptor<Md, detail::ORIENT_VEC_ROW> Dr;
    typedef detail::_ImplMtxAdaptor<Ms, detail::ORIENT_VEC_COL> Sc;
    typedef detail::_ImplMtxAdaptor<Ms, detail::ORIENT_VEC_ROW> Sr;
    
    // if a dimension match can be made with either orientation, proceed with the copy
    if ((Dc::rows(*into) == Sc::rows(src) and Dc::cols(*into) == Sc::cols(src)) or 
        (Dr::rows(*into) == Sr::rows(src) and Dr::cols(*into) == Sr::cols(src))) {
        detail::_mtxcopy(into, src);
    } else {
        throw DimensionMismatchException(Dc::rows(*into), Dc::cols(*into), Sc::rows(src), Sc::cols(src));
    }
#else
    detail::_mtxcopy(into, src);
#endif
}


/************************************
 * Matrix mul                       *
 ************************************/

/** 
 * matrix * matrix and matrix * vector multiplication
 * 
 * This function handles all multiplication between any two types of matrix,
 * or between any matrix and a vector. Thus, `Ma`, `Mb`, (the operands) and `Md` 
 * (the destination), may each be either a matrix or a vector.
 * If the left operand is a vector, it is assumed to be a row vector,
 * whereas right vector operands are column vectors.
 * 
 * This function ensures that all compile-time checkable dimensions
 * agree-- that is to say, that `(a x b) * (b x c) = (a x c)` holds.
 * If any object has a dynamic size, then the check will be deferred to
 * runtime, and a DimensionMismatchException thrown if the check fails. A 
 * compile-time dimension mismatch implies the program is invalid and compilation 
 * will error.
 * 
 * @param [out] into A writeable matrix or vector with dimensions `(a.rows() x b.cols())`
 * @param [in]  a    A matrix object or vector with dimension `(N x b.rows())`
 * @param [in]  b    A matrix object or vector with dimension `(a.cols() x N)`
 **/
#ifdef PARSING_DOXYGEN
template <typename Md, typename Ma, typename Mb> Md& mul(Md *into, const Ma &a, const Mb &b) {}
#endif

// (a x b) * (b x c) -> (a x c)
/* Implementation:
 * 
 * This function's main activity is to guarantee dimension safety and 
 * protect against memory aliasing. The actual calculation is performed
 * by _ImplMatrixMul, for which there are specializations over meaningful
 * operand types. _ImplMatrixMul also provides a resultant type (usually
 * a SimpleMatrix) for use when a temporary destination buffer is needed.
 * 
 * Because Vecs do not have Matrix-like 2D dimension, an _ImplMtxAdaptor
 * is used to make Vecs interoperable with matrices. The adaptor must be 
 * told whether a vector is expected to be a row vector or a column vector. 
 * Result vectors may be either a row or column vector depending on the 
 * order of the operands, so an _ImplVecMulOrient is used to select the 
 * appropriate orientation for result vectors.
 */
template <typename Md, typename Ma, typename Mb>
typename boost::enable_if_c<
            detail::MatrixMultipliable<Ma,Mb>::val and 
            detail::_ImplMtxResult<Ma, Mb, Md>::agreement == detail::MTX_RESULT_MATCH, 
        Md&>::type
mul(Md *into, const Ma &a, const Mb &b) {
    // adaptor types for each argument:
    typedef detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW> A;
    typedef detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL> B;
    typedef detail::_ImplMtxAdaptor<Md, detail::_ImplVecMulOrient<Ma,Mb,Md>::orient> D;
    // multiplier implementation:
    typedef detail::_ImplMtxMul<Ma,Mb> mult_t;
    typedef typename mult_t::return_t buffer_t;
    
    #ifdef GEOMC_MTX_CHECK_DIMS
    // do the source matrix dimensions agree?
    // any compiler worth its salt should eliminate this test for template instantiations with static-dimensioned operands
    if ((A::COLDIM == DYNAMIC_DIM or B::ROWDIM == DYNAMIC_DIM) and A::cols(a) != B::rows(b)) {
        throw DimensionMismatchException(A::rows(a), A::cols(a), B::rows(b), B::cols(b));
    }
    // does the destination matrix have correct dims?
    // again, runtime checks are to be avoided.
    if (((D::ROWDIM == DYNAMIC_DIM or A::ROWDIM == DYNAMIC_DIM) and A::rows(a) != D::rows(*into)) or 
        ((D::COLDIM == DYNAMIC_DIM or B::COLDIM == DYNAMIC_DIM) and B::cols(b) != D::cols(*into))) {
        throw DimensionMismatchException(D::rows(*into), D::cols(*into), A::rows(a), B::cols(b));
    }
    #endif
    
    #ifdef GEOMC_MTX_CHECK_ALIASING
    // allocate a temp destination matrix, in case of storage aliasing
    if (mtx_aliases_storage(*into, a) or mtx_aliases_storage(*into, b)) {
        buffer_t tmp = detail::_ImplMtxInstance<buffer_t>::instance(A::rows(a), B::cols(b));
        mult_t::mul(&tmp, a, b);
        detail::_mtxcopy(into, tmp);
    } else {
        mult_t::mul(into, a, b);
    }
    #else
    mult_t::mul(into, a, b);
    #endif
    
    return *into;
}

// vec <- mtx * mtx
// Special case where matrices have dynamic size and destination is a vector.
// The desired row/col orientation of the vector cannot be determined at compile
// time, so special logic is needed to check both configurations.
template <typename Md, typename Ma, typename Mb>
typename boost::enable_if_c<
                     detail::MatrixMultipliable<Ma,Mb>::val and 
                     detail::_ImplMtxResult<Ma,Mb,Md>::agreement == detail::MTX_RESULT_UNKNOWN,
            Md&>::type
mul(Md *into, const Ma &a, const Mb &b) {
    // adaptors types for arguments:
    typedef detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW> A;
    typedef detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL> B;
     //multiplier implementation
    typedef detail::_ImplMtxMul<Ma,Mb> mult_t;
    typedef typename mult_t::return_t buffer_t;
    
    #ifdef GEOMC_MTX_CHECK_DIMS
    // trial orientation adaptors:
    typedef detail::_ImplMtxAdaptor<Md, detail::ORIENT_VEC_ROW> Dr;
    typedef detail::_ImplMtxAdaptor<Md, detail::ORIENT_VEC_COL> Dc;
    
    // if neither configuration yields dimension agreement:
    if (not ((Dr::rows(*into) == A::rows(a) and Dr::cols(*into) == B::cols(b)) or
             (Dc::rows(*into) == A::rows(a) and Dc::cols(*into) == B::cols(b))) ) {
        throw DimensionMismatchException(Dr::rows(*into), Dr::cols(*into), A::rows(a), B::cols(b));
    }
    #endif
    
    #ifdef GEOMC_MTX_CHECK_ALIASING
    if (mtx_aliases_storage(*into, a) or mtx_aliases_storage(*into, b)) {
        buffer_t tmp = detail::_ImplMtxInstance<buffer_t>::instance(A::rows(a), B::cols(b));
        mult_t::mul(&tmp, a, b);
        detail::_mtxcopy(into, tmp);
    } else {
        mult_t::mul(into, a, b);
    }
    #else
    mult_t::mul(into, a, b);
    #endif
    
    return *into;
}

/**
 * Matrix multiplication.
 * 
 * `Ma` and `Mb` may either be matrices or vectors. Left operands are, if 
 * vectors, assumed to be rows, while right operands will be treated as columns.
 * The dimensions of `a` and `b` must satisfy `(a x b) * (b x c)`. If the dimensions
 * can be determined to mismatch at compile time, the program is considered invalid and 
 * will not compile. If either argument has dynamic size, the dimension check
 * will be performed at runtime, raising a `DimensionMismatchException` if it fails.
 * 
 * The return type will be chosen appropriately based on the arguments, and
 * will have dimension `(a x c)`.
 * 
 * @param [in]  a    A matrix object or vector with dimension `(N x b.rows())`
 * @param [in]  b    A matrix object or vector with dimension `(a.cols() x N)`
 * @return A new matrix or vector object containing the result of `a * b`. 
*/
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> Md mul(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
typename boost::enable_if_c<
            detail::MatrixMultipliable<Ma,Mb>::val, 
            typename detail::_ImplMtxMul<Ma,Mb>::return_t>::type
mul(const Ma &a, const Mb &b) {
    // adaptor types for operands:
    typedef detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW> A;
    typedef detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL> B;
    // other types:
    typedef detail::_ImplMtxMul<Ma,Mb> mult_t;
    typedef typename mult_t::return_t return_t;
    
    #ifdef GEOMC_MTX_CHECK_DIMS
    // dynamic dimensions match?
    if ((A::COLDIM == DYNAMIC_DIM or B::ROWDIM == DYNAMIC_DIM) and A::cols(a) != B::rows(b)) {
        throw DimensionMismatchException(A::rows(a), A::cols(a), B::rows(b), B::cols(b));
    }
    #endif
    
    return_t dest = detail::_ImplMtxInstance<return_t>::instance(A::rows(a), B::cols(b));
    mult_t::mul(&dest, a, b);
    return dest;
}

/*********************************
 * Matrix transpose              *
 *********************************/

/**
 * Matrix transpose.
 * 
 * @param [out] into A writeable matrix with dimensions `(m.cols(), m.rows())`
 * @param [in] m A matrix object
 */
#ifdef PARSING_DOXYGEN
template <typename Md, typename Mx> void transpose(Md *into, const Mx &m) {}
#endif
template <typename Md, typename Mx>
typename boost::enable_if_c<
    detail::IsMatrix<Md>::val and detail::IsMatrix<Mx>::val and
                  (Md::ROWDIM == Mx::COLDIM or
                   Md::ROWDIM *  Mx::COLDIM == 0) and
                  (Md::COLDIM == Mx::ROWDIM or 
                   Md::COLDIM *  Mx::ROWDIM == 0),
    void
>::type
transpose(Md *into, const Mx &m) {
    
    typedef detail::_ImplMtxTxpose<Mx> txpose_t;
    typedef typename txpose_t::return_t return_t;
    
    #ifdef GEOMC_MTX_CHECK_DIMS
    // dimension mismatch?
    if ((Md::ROWDIM * Mx::COLDIM == 0 and into->rows() != m.cols()) or 
        (Md::COLDIM * Mx::ROWDIM == 0 and into->cols() != m.rows())) {
        throw DimensionMismatchException(into->rows(), into->cols(), m.cols(), m.rows());
    }
    #endif
    
    #ifdef GEOMC_MTX_CHECK_ALIASING
    // memory aliasing?
    if (mtx_aliases_storage(*into, m)) {
        return_t buf = detail::_ImplMtxInstance<return_t>::instance(m.cols(), m.rows());
        txpose_t::transpose(&buf, m);
        detail::_mtxcopy(into, buf);
        return *into;
    }
    #endif
    txpose_t::transpose(into, m);
}

/**
 * Matrix transpose.
 * 
 * @param [in] m A matrix object.
 * @returns A transposed copy of `m`, of type appropriate for the argument, usually
 * a `SimpleMatrix`.
 */
#ifdef PARSING_DOXYGEN
template <typename Mx> Md transpose(const Mx &m) {}
#endif
template <typename Mx>
typename boost::enable_if_c<
            detail::IsMatrix<Mx>::val, 
            typename detail::_ImplMtxTxpose<Mx>::return_t>::type
transpose(const Mx &m) {
    typedef detail::_ImplMtxTxpose<Mx> txpose_t;
    typedef typename txpose_t::return_t return_t;
    
    return_t buf = detail::_ImplMtxInstance<return_t>::instance(m.cols(), m.rows());
    txpose_t::transpose(&buf, m);
    return buf;
}


/*********************************
 * Matrix inverse                *
 *********************************/

/**
 * Matrix inversion. `src` and `into` must be square matrices of the same dimension.
 * If a runtime check for square dimensions fails, a `NonsquareMatrixException` is raised.
 * 
 * @param [out] into A writeable matrix with dimensions equal to `src`.
 * @param [in]  src  A square matrix.
 * 
 * @return `false` if the matrix is singular and could not be inverted, `true` otherwise.
 */
#ifdef PARSING_DOXYGEN
template <typename Md, typename Mx> bool inv(Md *into, const Mx &src) {}
#endif
template <typename Md, typename Mx>
bool inv(Md *into, const Mx &src,
         M_ENABLE_IF_C(
             (detail::MatrixDimensionMatch<Md,Mx>::isStaticMatch and
             (Mx::ROWDIM == Mx::COLDIM or Mx::ROWDIM * Mx::COLDIM == DYNAMIC_DIM))
         )) {
    typedef detail::_ImplMtxInv<Mx> inv_t;
#ifdef GEOMC_MTX_CHECK_DIMS
    if (Mx::ROWDIM * Mx::COLDIM == DYNAMIC_DIM and src.rows() != src.cols()) {
        throw NonsquareMatrixException(src.rows(), src.cols());
    } 
    detail::MatrixDimensionMatch<Md,Mx>::check(*into, src);
#endif
    // matrix inv() implementations, unlike mul() or txpose(), shall assume
    // that the destination and source may be aliased. this design decision
    // was made because several inv() implementations are agnostic about
    // aliasing (they perform internal copies), making a check unnecessary.
    // therefore we delegate the aliasing check to those inv()s which care. 
    return inv_t::inv(into, src);
}


/**
 * Matrix inversion.
 * @param [in] m A square matrix.
 * @param [out] success Will be set to `false` if the matrix was singular and could not
 * be inverted, otherwise will be set to `true`.
 * @return A new matrix containing the inverse of `m`, or undefined data if 
 * `m` could not be inverted.
 */
#ifdef PARSING_DOXYGEN
template <typename Mx> Md inv(const Mx &m, bool *success) {}
#endif
template <typename Mx>
typename detail::_ImplMtxInv<Mx>::return_t inv(const Mx &m, bool *success, 
                M_ENABLE_IF_C(Mx::ROWDIM == Mx::COLDIM or Mx::ROWDIM * Mx::COLDIM == DYNAMIC_DIM)) {
    typedef detail::_ImplMtxInv<Mx> inv_t;
    typedef typename inv_t::return_t return_t;
    
#ifdef GEOMC_MTX_CHECK_DIMS
    if ((Mx::ROWDIM == DYNAMIC_DIM or Mx::COLDIM == DYNAMIC_DIM) and m.rows() != m.cols()) {
        throw NonsquareMatrixException(m.rows(), m.cols());
    }
#endif
    
    return_t into = detail::_ImplMtxInstance<return_t>::instance(m.rows(), m.cols());
    *success = inv_t::inv(&into, m);
    return into;
}


/*********************************
 * Matrix add/sub                *
 *********************************/

/**
 * Matrix addition. Add the corresponding elements of `a` and `b`, whose dimensions
 * must match.
 * 
 * @param [out] d A writeable matrix, whose dimensions must match `a` and `b`.
 * @param [in] a A matrix object
 * @param [in] b A matrix object
 */
#ifdef PARSING_DOXYGEN
template <typename Md, typename Ma, typename Mb> void add(Md *d, const Ma &a, const Mb &b) {}
#endif
template <typename Md, typename Ma, typename Mb>
typename boost::enable_if_c<
    detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch and
    detail::MatrixDimensionMatch<Ma,Md>::isStaticMatch,
    void
>::type 
add(Md *d, const Ma &a, const Mb &b) {
    typedef detail::_ImplMatrixAdd<Ma,Mb> add_t;
#ifdef GEOMC_MTX_CHECK_DIMS
    detail::MatrixDimensionMatch<Ma,Mb>::check(a,b);
    detail::MatrixDimensionMatch<Md,Ma>::check(*d,a);
#endif
    add_t::add(d,a,b);
}

/**
 * Matrix subtraction. Subtract the corresponding elements of `b` from `a`'s. The 
 * dimensions of `a` and `b` must match.
 * 
 * @param [out] d A writeable matrix, whose dimensions must match `a` and `b`.
 * @param [in] a A matrix object
 * @param [in] b A matrix object
 */
#ifdef PARSING_DOXYGEN
template <typename Md, typename Ma, typename Mb> void sub(Md *d, const Ma &a, const Mb &b) {}
#endif
template <typename Md, typename Ma, typename Mb>
typename boost::enable_if_c<
    detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch and
    detail::MatrixDimensionMatch<Ma,Md>::isStaticMatch,
    void
>::type 
sub(Md *d, const Ma &a, const Mb &b) {
    typedef detail::_ImplMatrixAdd<Ma,Mb> add_t;
#ifdef GEOMC_MTX_CHECK_DIMS
    detail::MatrixDimensionMatch<Ma,Mb>::check(a,b);
    detail::MatrixDimensionMatch<Md,Ma>::check(*d,a);
#endif
    add_t::sub(d,a,b);
}

/**
 * Matrix addition. Add the corresponding elements of `a` and `b`, whose
 * dimensions must match.
 * 
 * @param [in] a A matrix object
 * @param [in] b A matrix object
 * @return A new matrix containing `a + b`, usually a `SimpleMatrix`.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> Md add(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
typename boost::enable_if_c<
    detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch,
    typename detail::_ImplMatrixAdd<Ma,Mb>::return_t
>::type 
add(const Ma &a, const Mb &b) {
    typedef detail::_ImplMatrixAdd<Ma,Mb> add_t;
    typedef typename add_t::return_t return_t;
#ifdef GEOMC_MTX_CHECK_DIMS
    detail::MatrixDimensionMatch<Ma,Mb>::check(a,b);
#endif
    return_t into = detail::_ImplMtxInstance<return_t>::instance(a.rows(), a.cols());
    add_t::add(&into, a, b);
    return into;
}

/**
 * Matrix subtraction. Subtract the corresponding elements of `b` from `a`'s. The
 * dimensions of `a` and `b` must match.
 * 
 * @param [in] a Matrix object
 * @param [in] b Matrix object
 * @return A new matrix containing `a - b`, usually a `SimpleMatrix`.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> Md sub(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
typename boost::enable_if_c<
    detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch,
    typename detail::_ImplMatrixAdd<Ma,Mb>::return_t
>::type 
sub(const Ma &a, const Mb &b) {
    typedef detail::_ImplMatrixAdd<Ma,Mb> add_t;
    typedef typename add_t::return_t return_t;
#ifdef GEOMC_MTX_CHECK_DIMS
    detail::MatrixDimensionMatch<Ma,Mb>::check(a,b);
#endif
    return_t into = return_t(a.rows(), a.cols());
    add_t::sub(&into, a, b);
    return into;
}


/*********************************
 * Matrix scale                  *
 *********************************/

/**
 * Scalar muliplication on matrices. In other words, multiply all the elements
 * of `m` by scalar value `k`.
 * 
 * @param [out] d A writeable matrix, whose dimensions must match those of `m`.
 * @param [in]  k Scalar constant (whose type satisfies `boost::is_scalar<U>`).
 * @param [in]  m Matrix object to be scaled.
 */
#ifdef PARSING_DOXYGEN
template <typename U typename Mx, typename Md> void scale(Md *d, U k, const Mx &m) {}
#endif
template <typename U, typename Mx, typename Md>
typename boost::enable_if_c<
    detail::IsMatrix<Mx>::val and
    boost::is_scalar<U>::value and
    detail::MatrixDimensionMatch<Mx,Md>::isStaticMatch,
    void
>::type 
scale(Md *d, U k, const Mx &m) {
    typedef detail::_ImplMatrixScale<Mx> scale_t;
    typedef typename scale_t::return_t return_t;
#ifdef GEOMC_MTX_CHECK_DIMS
    detail::MatrixDimensionMatch<Mx,Md>::check(m,*d);
#endif
    scale_t::scale(d, k, m);
}

/**
 * Scalar muliplication on matrices. In other words, multiply all the elements
 * of `m` by scalar value `k`.
 * 
 * @param [in]  k Scalar constant (whose type satisfies `boost::is_scalar<U>`).
 * @param [in]  m Matrix object to be scaled.
 * @return A scaled copy of `m`.
 */
#ifdef PARSING_DOXYGEN
template <typename U typename Mx> Mx scale(U k, const Mx &m) {}
#endif
template <typename U, typename Mx>
typename boost::enable_if_c<
    detail::IsMatrix<Mx>::val and
    boost::is_scalar<U>::value,
    typename detail::_ImplMatrixScale<Mx>::return_t
>::type 
scale(U k, const Mx &m) {
    typedef detail::_ImplMatrixScale<Mx> scale_t;
    typedef typename scale_t::return_t return_t;
    return_t into = return_t(m.rows(), m.cols());
    scale_t::scale(&into, k, m);
    return into;
}

/*********************************
 * Matrix operators              *
 *********************************/

/**
 * scalar * matrix
 * 
 * All elements of `m` are multiplied by `k`. `Mx` must be a matrix type, and
 * `U` must satisfy `boost::is_scalar<U>`.
 */
#ifdef PARSING_DOXYGEN
template <typename U, typename Mx> Mx operator*(U k, const Mx &m) {}
#endif
template <typename U, typename Mx>
inline typename boost::enable_if_c<
    detail::IsMatrix<Mx>::val and boost::is_scalar<U>::value,
    typename detail::_ImplMatrixScale<Mx>::return_t
>::type operator*(U k, const Mx &m) {
    return scale(k,m);
}

/**
 * matrix * scalar
 * 
 * All elements of `m` are multiplied by `k`. `Mx` must be a matrix type, and
 * `U` must satisfy `boost::is_scalar<U>`.
 */
#ifdef PARSING_DOXYGEN
 template <typename U, typename Mx> Mx operator*(const Mx &m, U k) {}
#endif
template <typename U, typename Mx>
inline typename boost::enable_if_c<
    detail::IsMatrix<Mx>::val and boost::is_scalar<U>::value,
    typename detail::_ImplMatrixScale<Mx>::return_t
>::type operator*(const Mx &m, U k) {
    return scale(k,m);
}

/**
 * matrix + matrix
 * 
 * Add the elements of `a` and `b`. The return type will be chosen appropriately
 * based on the arguments (usually a `SimpleMatrix`). `Ma` and `Mb` must both be 
 * matrix types.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> Mx operator+(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
inline typename detail::_ImplMatrixAddReturnType<Ma,Mb>::return_t
operator+(const Ma &a, const Mb &b) {
    return add(a,b);
}

/**
 * matrix - matrix. 
 * 
 * Subtract the elements of `b` from `a`. The return type will be chosen appropriately
 * based on the arguments (usually a `SimpleMatrix`). `Ma` and `Mb` must both be
 * matrix types.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> Mx operator-(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
inline typename detail::_ImplMatrixAddReturnType<Ma,Mb>::return_t
operator-(const Ma &a, const Mb &b) {
    return sub(a,b);
}

/**
 * matrix * matrix
 * 
 * Matrix multiplication is performed on `a` and `b`.
 * 
 * `Ma` and `Mb` may either be matrices or vectors (but may not _both_ be vectors;
 * this is handled by a different multiplication operator). Left operands are, if 
 * vectors, assumed to be rows, while right operands will be treated as columns.
 * The dimensions of `a` and `b` must satisfy `(a x b) * (b x c)`. If the dimensions
 * can be determined to mismatch at compile time, the program is considered invalid and 
 * will not compile. If either argument has dynamic size, the dimension check
 * will be performed at runtime, raising a `DimensionMismatchException` if it fails.
 * 
 * The return type will be chosen appropriately based on the arguments, and
 * will have dimension `(a x c)`.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> Mx operator*(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
inline typename detail::MatrixMultReturnType<Ma,Mb>::return_t
operator*(const Ma &a, const Mb &b) {
    return mul<Ma,Mb>(a, b);
}

/**
 * matrix == matrix. 
 * 
 * Matrices `a` and `b` are equal if and only if `a` and `b` have the same dimension 
 * and all corresponding elements are equal. `Ma` and `Mb` must both be matrix types. 
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> bool operator==(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
typename boost::enable_if_c<
        detail::IsMatrix<Ma>::val and detail::IsMatrix<Mb>::val and
        detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch, 
    bool>::type
operator==(const Ma &a, const Mb &b) {
    if (not detail::MatrixDimensionMatch<Ma,Mb>::isMatch(a,b)) {
        // dimension mismatch, cannot be the same.
        return false;
    }
    
    return detail::mtxequal(a,b);
}

// mtx == mtx 
// (static dimension mismatch; never equal)
template <typename Ma, typename Mb>
typename boost::enable_if_c<
                detail::IsMatrix<Ma>::val and detail::IsMatrix<Mb>::val and
            not detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch, 
        bool>::type
operator==(const Ma &a, const Mb &b) {
    return false;
}

/**
 * matrix != matrix. 
 * 
 * Matrices `a` and `b` are unequal unless `a` and `b` have the same dimension 
 * and all corresponding elements are equal.
 */
#ifdef PARSING_DOXYGEN
template <typename Ma, typename Mb> bool operator!=(const Ma &a, const Mb &b) {}
#endif
template <typename Ma, typename Mb>
inline typename boost::enable_if_c<detail::IsMatrix<Ma>::val and detail::IsMatrix<Mb>::val, bool>::type
operator!=(const Ma &a, const Mb &b) {
    return not (a == b);
}

#ifdef GEOMC_LINALG_USE_STREAMS

/**
 * stream << matrix 
 * 
 * Output a complete textual representation of `m` to stream `s`. Elements are
 * separated by spaces and rows are separated by newlines.
 */
#ifdef PARSING_DOXYGEN
template <typename Mx> std::ostream& operator<< (std::ostream &s, const Mx &m) {}
#endif
template <typename Mx>
typename boost::enable_if_c<detail::IsMatrix<Mx>::val, std::ostream &>::type
operator<<(std::ostream &s, const Mx &mtx) {
    s << std::setfill(' '); //xxx statefulness bad
    s << "[ ";
    for (int r = 0; r < mtx.rows(); r++){
        if (r != 0) s << "  ";
        for (int c = 0; c < mtx.cols(); c++){
            s << std::setw(9) << mtx.get(r,c) << " ";
        }
        if (r == mtx.rows() - 1){
            // last row, close the matrix.
            s << "]";
        }
        s << std::endl;
    }
    return s;
}

#endif

}; // namespace geom

/// @} //ingroup linalg

#endif /* MATRIX_H_ */
