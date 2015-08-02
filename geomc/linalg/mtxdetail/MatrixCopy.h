/* 
 * File:   MatrixCopy.h
 * Author: tbabb
 *
 * Created on September 14, 2013, 6:56 PM
 */

#ifndef MATRIXCOPY_H
#define	MATRIXCOPY_H

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>


namespace geom {
namespace detail {

    template <typename Md, typename Mx>
    inline void _mtxcopy(Md *into, const Mx &src) {
        std::copy(src.begin(), src.end(), into->begin()); //also works for vectors. neat!
    }


    template <index_t N>
    inline void _mtxcopy(geom::PermutationMatrix<N> *into, const geom::PermutationMatrix<N> &src) {
        into->setRowSources(src.getRowSources());
    }


    template <typename Td, index_t Md, index_t Nd, 
              typename Ts, index_t Ms, index_t Ns>
    void _mtxcopy(geom::DiagMatrix<Td, Md, Nd> *into, const geom::DiagMatrix<Ts, Ms, Ns> &src) {
        std::copy(src.diagonal_begin(), src.diagonal_end(), into->diagonal_begin());
    }


}; // namespace detail


/**
 * @addtogroup matrix
 * @{
 */

/*********************************
 * User-facing functions         *
 *********************************/

/** 
 * Copy the contents of `src` into `into`. This function may be optimized to 
 * perform better than running `std::copy()` on some matrix types' iterators. It
 * will also succeed in certain situations where `std::copy()` would fail (for
 * example when copying the contents of a diagonal matrix to another, any attempt
 * to write off the diagonal will throw an exception).
 * 
 * `Matrix` and `Matrix1` must be matrix or vector types whose dimensions match.
 * If the dimensions can be determined to mismatch at compile-time, the program
 * is considered invalid and the compilation will fail. If either object has dynamic
 * size, the check will be deferred to runtime, throwing a `DimensionMismatchException`
 * if the check fails.
 * 
 * @param [out] into A writeable matrix with dimensions matching `src`'s.
 * @param [in]  src  A matrix object.
 */
#ifdef PARSING_DOXYGEN
template <typename Matrix, typename Matrix1> void mtxcopy(Matrix *into, const Matrix1 &src) {}
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
    if ((D::ROWDIM * S::ROWDIM == DYNAMIC_DIM and D::rows(*into) != S::rows(src)) or
        (D::COLDIM * S::COLDIM == DYNAMIC_DIM and D::cols(*into) != S::cols(src))) {
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

/// @} // addtogroup matrix

}; // namespace geom

#endif	/* MATRIXCOPY_H */

