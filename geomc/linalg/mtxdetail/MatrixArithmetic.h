/* 
 * File:   MatrixArithmetic.h
 * Author: tbabb
 *
 * Created on July 28, 2013, 11:53 AM
 */

#ifndef MATRIXARITHMETIC_H
#define	MATRIXARITHMETIC_H

#include <boost/utility/enable_if.hpp>
#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>

namespace geom {
namespace detail {
    
/************************************
 * Matrix add/sub                   *
 ************************************/

template <typename Ma, typename Mb>
struct _ImplMatrixAdd {
    
    //todo: should use decltype(Ma::elem_t + Mb::elem_t)
    typedef SimpleMatrix<typename Ma::elem_t, 
                         Ma::ROWDIM * Mb::COLDIM == 0 ? 0 : Ma::ROWDIM, 
                         Ma::COLDIM * Mb::COLDIM == 0 ? 0 : Ma::COLDIM>
            return_t;
    
    template <typename Md>
    static void add(Md *d, const Ma &a, const Mb &b) {
        typename Md::iterator d_i = d->begin();
        typename Ma::const_iterator a_i = a.begin();
        typename Mb::const_iterator b_i = b.begin();
        for (; d_i != d->end(); ++d_i, ++a_i, ++b_i) {
            *d_i = (*a_i) + (*b_i);
        }
    }
    
    template <typename Md>
    static void sub(Md *d, const Ma &a, const Mb &b) {
        typename Md::iterator d_i = d->begin();
        typename Ma::const_iterator a_i = a.begin();
        typename Mb::const_iterator b_i = b.begin();
        for (; d_i != d->end(); ++d_i, ++a_i, ++b_i) {
            *d_i = (*a_i) - (*b_i);
        }
    }
    
};

//////////// Diagonal case ////////////

template <typename T, index_t M, index_t N,
          typename S, index_t J, index_t K>
struct _ImplMatrixAdd < DiagMatrix<T,M,N>, DiagMatrix<S,J,K> > {
    
    typedef DiagMatrix<T,M,N> Ma;
    typedef DiagMatrix<S,J,K> Mb;
    
    //TODO: use decltype(T+S)
    typedef DiagMatrix<T, 
                       M * J == 0 ? 0 : M,
                       N * K == 0 ? 0 : N>
            return_t;
    
    // arbitrary destination
    
    template <typename Md>
    static void add(Md *dest, const Ma &a, const Mb &b) {
        dest->setZero();
        const index_t diag = std::min(a.rows(), a.cols());
        
        const T *a_i = a.diagonal_begin();
        const S *b_i = b.diagonal_begin();
        for (index_t i = 0; i < diag; ++a_i, ++b_i, ++i) {
            dest->set(i, i, *a_i + *b_i);
        }
    }
    
    template <typename Md>
    static void sub(Md *dest, const Ma &a, const Mb &b) {
        dest->setZero();
        const index_t diag = std::min(a.rows(), a.cols());
        
        const T *a_i = a.diagonal_begin();
        const S *b_i = b.diagonal_begin();
        for (index_t i = 0; i < diag; ++a_i, ++b_i, ++i) {
            dest->set(i,i, *a_i - *b_i);
        }
    }
    
    // diagonal destination
    // (save a setZero, and use bare pointers).
    
    template <typename U, index_t P, index_t Q>
    static void add(DiagMatrix<U,P,Q> *dest, const Ma &a, const Mb &b) {
        const T *a_i = a.diagonal_begin();
        const S *b_i = b.diagonal_begin();
        U *d_i = dest->diagonal_begin();
        
        for (; a_i != a.diagonal_end(); ++a_i, ++b_i, ++d_i) {
            *d_i = *a_i + *b_i;
        }
    }
    
    template <typename U, index_t P, index_t Q>
    static void sub(DiagMatrix<U,P,Q> *dest, const Ma &a, const Mb &b) {
        const T *a_i = a.diagonal_begin();
        const S *b_i = b.diagonal_begin();
        U *d_i = dest->diagonal_begin();
        
        for (; a_i != a.diagonal_end(); ++a_i, ++b_i, ++d_i) {
            *d_i = *a_i - *b_i;
        }
    }
    
};

    
// stupid stupid workaround because c++ is a stupid stupid language:
template <typename Ma, typename Mb, typename Enable=void>
struct _ImplMatrixAddReturnType {
    // no return type. engage SFINAE.
};

// when using enable_if<bool, _implmtxadd<a,b>::return_t> elsewhere,
// compilation will fail because even if <bool> is false and "protects" 
// against instantiation, c++ will try to expand the argument type. if a or b 
// are non-matrices, then the expansion will fail. This is not SFINAE because by 
// the time we are evaluating what ::return_t is, we are no longer substituting,
// we are type-expanding. I believe this is a flaw in the c++ language design.

template <typename Ma, typename Mb>
struct _ImplMatrixAddReturnType <Ma, Mb, 
            typename boost::enable_if_c<
                detail::IsMatrix<Ma>::val and
                detail::IsMatrix<Mb>::val and
                detail::MatrixDimensionMatch<Ma,Mb>::isStaticMatch,
                void
            >::type > {
    typedef typename _ImplMatrixAdd<Ma,Mb>::return_t return_t;
};
    
/************************************
 * Matrix scalar mul                *
 ************************************/

template <typename Mx>
struct _ImplMatrixScale {
    typedef Mx return_t;
    
    template <typename Md, typename U>
    static void scale(Md *d, U k, const Mx &m) {
        for (index_t r = 0; r < m.rows(); ++r) {
            for (index_t c = 0; c < m.cols(); ++c) {
                d->set(r, c, k * m.get(r,c));
            }
        }
    }
};

////////// diagonal case //////////

template <typename T, index_t M, index_t N>
struct _ImplMatrixScale < DiagMatrix<T,M,N> > {
    typedef DiagMatrix<T,M,N> return_t;
    
    template <typename Md, typename U>
    static void scale(Md *d, U k, const DiagMatrix<T,M,N> &m) {
        d->setZero();
        const index_t diag = std::min(m.rows(), m.cols());
        const T *m_i = m.diagonal_begin();
        for (index_t i = 0; i < diag; ++i, ++m_i) {
            d->set(i,i, k * (*m_i));
        }
    }
    
    template <typename S, index_t J, index_t K, typename U>
    static void scale(DiagMatrix<S,J,K> *d, U k, const DiagMatrix<T,M,N> &m) {
        S *d_i = d->diagonal_begin();
        const T *m_i = m.diagonal_begin();
        for (; m_i != m.diagonal_end(); ++d_i, ++m_i) {
            *d_i = k * (*m_i);
        }
    }
};


};  // end namespace detail
};  // end namespace geom


#endif	/* MATRIXARITHMETIC_H */

