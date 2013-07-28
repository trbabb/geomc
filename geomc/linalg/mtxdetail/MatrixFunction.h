/*
 * MatrixFunction.h
 *
 *  Created on: May 23, 2013
 *      Author: tbabb
 */


#ifndef MATRIXFUNCTION_H_
#define MATRIXFUNCTION_H_

#include <boost/type_traits/is_scalar.hpp>
#include <geomc/linalg/mtxdetail/MatrixFunctionImpl.h>
#include <geomc/linalg/mtxdetail/MatrixMult.h>
#include <geomc/linalg/mtxdetail/LUDecomp.h>
#include <geomc/linalg/mtxdetail/MatrixInv.h>
#include <geomc/linalg/mtxdetail/MatrixArithmetic.h>

#ifdef GEOMC_LINALG_USE_STREAMS
  #include <iostream>
  #include <iomanip>
#endif

namespace geom {

/************************************
 * Storage aliasing check           *
 ************************************/

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

/* Some combinations of matrices (like sparse <- permute, or
 * permute <- permute, or diag <- diag) can have copy
 * operations that are simpler/faster than element-wise copying.
 * This function abstracts the algorithm selection.
 */

template <typename Md, typename Ms>
void mtxcopy(Md *into, const Ms &src,
                typename boost::enable_if_c<
                     LINALG_DIM_AGREE(Md,Ms) and
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

/* Matrix <-> Matrix / Matrix <-> Vector mult operation
 * 
 * This function handles all multiplication between any two types of matrix,
 * or between any matrix and a vector. Thus, Ma, Mb, (the operands) and Md 
 * (the destination), may each be either a matrix or a vector.
 * If the left operand is a vector, it is assumed to be a row vector,
 * whereas right vector operands are column vectors.
 * 
 * This function ensures that all compile-time checkable dimensions
 * agree-- that is to say, that (a x b) * (b x c) -> (a x c) holds.
 * If any object has a dynamic size, then the check will be deferred to
 * runtime, and a DimensionMismatchException thrown if the check fails.
 * 
 * Type verification also guarantees that this function will not hide
 * other templated mul() functions; i.e. if the operands are not matrices
 * or vectors, this template will not match.
 * 
 * Implementation:
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
// mtx <- mtx * mtx
// (a x b) * (b x c) -> (a x c)
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

// matrix * matrix -> matrix
// return type chosen appropriately for arguments (e.g. diag * diag -> diag)    
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

// matrix <- matrix transpose
template <typename Md, typename Mx>
Md& transpose(Md *into, const Mx &m, 
              M_ENABLE_IF_C(
                  detail::IsMatrix<Md>::val and detail::IsMatrix<Mx>::val and
                  (Md::ROWDIM == Mx::COLDIM or
                   Md::ROWDIM *  Mx::COLDIM == 0) and
                  (Md::COLDIM == Mx::ROWDIM or 
                   Md::COLDIM *  Mx::ROWDIM == 0) )) {
    
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
    detail::_ImplMtxTxpose<Mx>(into, m);
    return *into;
}


// matrix transpose -> matrix
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

// mtx <- inv(mtx)
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


// inv(mtx) -> mtx
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

// mtx <- mtx + mtx
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

// mtx <- mtx - mtx
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

// mtx + mtx -> mtx
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

// mtx - mtx -> mtx
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

// mtx <- const * mtx
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

// const * mtx -> mtx
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

// constant * mtx
template <typename U, typename Mx>
inline typename boost::enable_if_c<
    detail::IsMatrix<Mx>::val and boost::is_scalar<U>::value,
    typename detail::_ImplMatrixScale<Mx>::return_t
>::type operator*(U k, const Mx &m) {
    return scale(k,m);
}

// mtx * constant
template <typename U, typename Mx>
inline typename boost::enable_if_c<
    detail::IsMatrix<Mx>::val and boost::is_scalar<U>::value,
    typename detail::_ImplMatrixScale<Mx>::return_t
>::type operator*(const Mx &m, U k) {
    return scale(k,m);
}

// mtx + mtx
template <typename Ma, typename Mb>
inline typename detail::_ImplMatrixAddReturnType<Ma,Mb>::return_t
operator+(const Ma &a, const Mb &b) {
    return add(a,b);
}

// mtx - mtx
template <typename Ma, typename Mb>
inline typename detail::_ImplMatrixAddReturnType<Ma,Mb>::return_t
operator-(const Ma &a, const Mb &b) {
    return sub(a,b);
}

// mtx * mtx
template <typename Ma, typename Mb>
//error if uncommented. cannot figure out why.
//typename boost::enable_if_c<(detail::MatrixMultipliable<Ma,Mb>::val), typename detail::_ImplMtxMul<Ma,Mb>::return_t>::type 
inline typename boost::enable_if_c<MATRIX_MUL_DIM_AGREE(Ma,Mb), typename detail::_ImplMtxMul<Ma,Mb>::return_t>::type 
operator*(const Ma &a, const Mb &b) {
    return mul<Ma,Mb>(a, b);
}

// mtx == mtx
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

// mtx != mtx
template <typename Ma, typename Mb>
inline typename boost::enable_if_c<detail::IsMatrix<Ma>::val and detail::IsMatrix<Mb>::val, bool>::type
operator!=(const Ma &a, const Mb &b) {
    return not (a == b);
}

#ifdef GEOMC_LINALG_USE_STREAMS
// stream output operator
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

#endif /* MATRIXFUNCTION_H_ */
