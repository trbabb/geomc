/*
 * LUDecomp.h
 *
 *  Created on: Jul 11, 2013
 *      Author: tbabb
 */

#ifndef LUDECOMP_H_
#define LUDECOMP_H_

#include "linalg/Matrix.h"

namespace geom {

namespace detail {

/************************************************
 * Matrix LU decomposition                      *
 *                                              *
 * Implementation using Doolittle's algorithm   *
 ************************************************/

#define _MxElem(r,c) m[cols*r + c]

// we operate on a bare data array (in row-contiguous order).
// this saves on multiple instantiations of this method for
// different static sizes of SimpleMatrix. It's 7% faster too!
template <typename T>
bool _ImplDecompPLU(T* m, index_t rows, index_t cols, index_t *reorder, bool *swap_parity) {
    const index_t n = std::min(rows, cols);
    *swap_parity = false;
    // fill permutation array
    for (index_t i = 0; i < n; i++){
        reorder[i] = i;
    }
    
    for (index_t i = 0; i < n - 1; i++) {
        // find pivot
        T biggest = std::abs(_MxElem(i,i));
        index_t pvt = i;
        for (index_t p = i + 1; p < rows; p++) {
            T pvt_val = std::abs(_MxElem(p,i));
            if (pvt_val > biggest) {
                pvt = p;
                biggest = pvt_val;
            }
        }
        
        if (biggest == 0) {
            // singular matrix
            // could test against an epsilon to 
            // find ill-conditioned matrices
            return false;
        } else if (pvt != i) {
            // swap row <i> with row <pvt>
            for (index_t c = 0; c < cols; c++) {
                std::swap(_MxElem(i,   c),
                          _MxElem(pvt, c));
            }
            // make a note of the permutation in <P>
            std::swap(reorder[i], reorder[pvt]);
            *swap_parity = !*swap_parity;
        }
        
        // eliminate lower elements
        T a = _MxElem(i,i);
        for (index_t r = i + 1; r < rows; r++) {
            T b = _MxElem(r,i) / a;
            for (index_t c = i + 1; c < cols; c++) {
                // R_r = R_r - b * R_i
                T src_elem = _MxElem(i,c);
                T dst_elem = _MxElem(r,c);
                _MxElem(r,c) = dst_elem - b * src_elem;
            }
            // set the lower matrix
            _MxElem(r,i) = b;
        }
    }
    
    return true;
}

#undef _MxElem

}; // namespace detail

//////////// PLU class ////////////

// LU is stores both the upper and lower triangular parts of the decomposition.
// The diagonal one elements (which belong to L) are not stored.

// P has dimension (U.rows() x U.rows())
// LU has the dimension of the source

template <typename T, index_t M, index_t N>
class PLUDecomposition {
    
protected:
    SimpleMatrix<T,M,N> LU;
    PermutationMatrix<M> P;
    bool singular;
    bool swap_parity;
    
    PLUDecomposition(index_t n_r, index_t n_c):
            LU(n_r, n_c),
            P(n_r),
            singular(false),
            swap_parity(false) {}
    
public:
    template <typename Mx>
    explicit PLUDecomposition(const Mx& m, 
                              typename boost::enable_if_c<detail::MatrixDimensionMatch<SimpleMatrix<T,M,N>, Mx>::isStaticMatch, int>::type dummy=0):
            LU(m.rows(), m.cols()),
            P(m.rows()),
            singular(false),
            swap_parity(false) {
        
        mtxcopy(&LU, m);
        // TODO: this alloc isn't necessary if we become a friend of PermutationMatrix
        detail::TemporaryStorage<index_t, Mx::ROWDIM> reorder(m.rows());
        bool ok = detail::_ImplDecompPLU(LU.begin(), m.rows(), m.cols(), reorder.get(), &swap_parity);

        if (not ok) {
            singular = true;
        } else {
            P.setRowSources(reorder.get());
        }
    }
    
public:
    
    const PermutationMatrix<M>& getP() const {
        return P;
    }
    
    const SimpleMatrix<T,M,N>& getLU() const {
        return LU;
    }
    
    const SimpleMatrix<T,M,M> getL() const {
        typedef detail::_ImplMtxInstance< SimpleMatrix<T,M,M> > instancer;
        SimpleMatrix<T,M,M> out = instancer::instance(LU.rows(), LU.rows());
        _copyL(out);
        return out;
    }
    
    const SimpleMatrix<T,M,N> getU() const {
        typedef detail::_ImplMtxInstance< SimpleMatrix<T,M,N> > instancer;
        SimpleMatrix<T,M,N> out = instancer::instance(LU.rows(), LU.cols());
        _copyU(out);
        return out;
    }
    
    inline void getL(SimpleMatrix<T,M,N> *into) const {
        _copyL(into);
    }
    
    inline void getU(SimpleMatrix<T,M,N> *into) const {
        _copyU(into);
    }
    
    // only vectors of dimension <U.rows()> are permitted.
    // we'll assume (perhaps unfairly?) that the client
    // is passing the right thing, since there is no
    // better way to enforce this.
    template <typename S>
    inline void linearSolve(S *dest, const S *b) const {
        
        #ifdef GEOMC_MTX_CHECK_DIMS
        _checkIsSquare();
        #endif
        
        #ifdef GEOMC_MTX_CHECK_ALIASING
        if (dest == b) {
            // because of the permutation, <b> will be destructively
            // updated as it is read.
            index_t n = LU.rows();
            detail::TemporaryStorage<S,M> buf(n);
            std::copy(b, b+n, buf.get());
            _linearSolve(dest, buf.get());
            return;
        }
        #endif
        
        _linearSolve(dest, b);
     }
    
    
    template <typename S, index_t K>
    inline Vec<S,K> linearSolve(const Vec<S,K> &b) const {
        Vec<S,K> dest;

        // vectors get a dimension check, because we can:
        #ifdef GEOMC_MTX_CHECK_DIMS
        if ((M == DYNAMIC_DIM or M != K) && LU.rows() != K) {
            throw DimensionMismatchException(LU.rows(), 1, K, 1);
        }
        _checkIsSquare();
        #endif
        
        _linearSolve(dest.begin(), b.begin());
        return dest;
    }
    
    
    template <typename S, index_t J, index_t K>
    void inverse(SimpleMatrix<S,J,K> *into) const {
        
        #ifdef GEOMC_MTX_CHECK_DIMS
        _checkIsSquare();
        // destination is valid?
        if ((J*K == 0 or J != M or K != N) and 
            (into->rows() != LU.rows() or into->cols() != LU.cols())) {
            throw DimensionMismatchException(into->rows(), into->cols(), LU.rows(), LU.cols());
        }
        #endif
        _inverseTranspose(into->begin());
        into->transpose();
    }
    
    
    T det() const {
        const index_t n = LU.rows();
        T k = getParity();
        for (index_t i = 0; i < n; i++) {
            k *= LU.get(i,i);
        }
        return k;
    }
    
    
    inline bool isSingular() const {
        return singular;
    }
    
    
    inline int getParity() const {
        return swap_parity ? -1 : 1;
    }
    
protected:
    
    template <typename Mx>
    void _copyL(Mx *into) const {
        for (index_t r = 0; r < LU.rows(); r++) {
            for (index_t c = 0; c < LU.cols(); c++) {
                T v;
                if (c < r) {
                    v = LU.get(r,c);
                } else if (c == r) {
                    v = 1;
                } else {
                    v = 0;
                }
                into->set(r,c,v);
            }
        }
    }
    
    template <typename Mx>
    void _copyU(Mx *into) const {
        for (index_t r = 0; r < LU.rows(); r++) {
            for (index_t c = 0; c < LU.cols(); c++) {
                T v;
                if (c >= r) {
                    v = LU.get(r,c);
                } else {
                    v = 0;
                }
                into->set(r,c,v);
            }
        }
    }
    
    inline void _checkIsSquare() const {
        if ((M * N == 0 or M != N) and LU.rows() != LU.cols()) {
            throw NonsquareMatrixException(LU.rows(), LU.cols());
        }
    }
    
    
    template <typename S>
    void _linearSolve(S *dest, const S *b, bool permute=true) const {
        const index_t  n = LU.rows();
        const index_t *p = P.getRowSources(); // permutation map
        
        // <y> and <dest>'s elements are used such that
        // y[i] is never read after dest[i] is written.
        // thus we may collapse their storage and save space:
        S *y = dest;

        // LUx = Pb
        // (Ux) is a vector, so let's solve for it:
        // Ly  = Pb
        y[0] = permute ? b[p[0]] : b[0];
        for (index_t r = 1; r < n; r++) {
            y[r] = permute ? b[p[r]] : b[r];
            for (index_t c = 0; c < r; c++) {
                y[r] -= y[c] * LU.get(r,c);
            }
        }
        // now with y, we may obtain x from:
        // Ux = y
        for (index_t r = n - 1; r >= 0; r--) {
            // dest[r] = y[r]; // nop; dest and y are the same!
            for (index_t c = n - 1; c > r; c--) {
                dest[r] -= dest[c] * LU.get(r,c);
            }
            dest[r] /= LU.get(r,r);
        }
    }
    
    
    template <typename S>
    void _inverseTranspose(S *dest) const {
        // Set LUx = PI and solve for x, choosing columns of PI one at a time.
        // Here, we use the rows of our destination matrix as though they are
        // column vectors, and the caller will transpose.
        const index_t *p = P.getColSources();
        const index_t n = LU.rows();
        std::fill(dest, dest + (n*n), 0);
        for (index_t i = 0; i < n; i++, dest += n) {
            dest[p[i]] = 1;
            _linearSolve(dest, dest, false);
        }
    }
};

};


#endif /* LUDECOMP_H_ */
