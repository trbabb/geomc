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

//TODO: do this only on bare data?
//      yes. all matrices get copied into one working copy, where the 
//      decomposition happens.
//TODO: destructive update by putting both L and U into
//      src matrix saves space and probably improves
//      cache performance/branch prediction
//TODO: could push storage down into PLU, and manipulate the "collapsed" matrix.
//      then add .getP(), .getL(), .getU()?
//      this would have the added benefit of reducing bloat
//      storage would be of dim M*N, and the P-matrix would remain as-is.

template <typename Ml, typename Mu>
bool _ImplDecompPLU(Ml *mat_lower, Mu *mat_upper_src, index_t *reorder, bool *swap_parity) {
    typedef typename Mu::elem_t T;
    typedef typename Ml::elem_t S;
    
    const index_t n_r = mat_upper_src->rows();
    const index_t n_c = mat_upper_src->cols();
    const index_t n = std::min(n_r, n_c);
    *swap_parity = false;
    // fill permutation array
    for (index_t i = 0; i < n; i++){
        reorder[i] = i;
    }
    
    for (index_t i = 0; i < n - 1; i++) {
        // find pivot
        T biggest = std::abs(mat_upper_src->get(i,i));
        index_t pvt = i;
        for (index_t p = i + 1; p < n_r; p++) {
            T pvt_val = std::abs(mat_upper_src->get(p,i));
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
            // everything to the left of col i is zero for these rows.
            for (index_t c = i; c < n_c; c++) {
                std::swap(mat_upper_src->get(i,   c),
                          mat_upper_src->get(pvt, c));
            }
            // ...and the lower matrix
            // everything at and beyond the pivot column remains in place
            for (index_t c = 0; c < i; c++) {
                std::swap(mat_lower->get(i,   c),
                          mat_lower->get(pvt, c));
            }
            // make a note of the permutation in <P>
            std::swap(reorder[i], reorder[pvt]);
            *swap_parity = !*swap_parity;
        }
        
        // eliminate lower elements
        T a = mat_upper_src->get(i,i);
        for (index_t r = i + 1; r < n_r; r++) {
            T b = mat_upper_src->get(r,i) / a;
            // targeted element goes to zero.
            mat_upper_src->set(r,i,0);
            for (index_t c = i + 1; c < n_c; c++) {
                // R_r = R_r - b * R_i
                T src_elem = mat_upper_src->get(i,c);
                T dst_elem = mat_upper_src->get(r,c);
                mat_upper_src->set(r, c, dst_elem - b * src_elem);
            }
            // set the lower matrix
            mat_lower->set(r, i, b);
        }
    }
    
    return true;
}

//////////// PLU class, base template ////////////

// base class necessary because matrix constructors have mandatory 
// arguments iff one or more dimensions are dynamic. this mandates
// special constructors for each case, and thus entire class templates
// for each case. did I mention that c++ is dumb?

// L and P always have dimension (U.rows() x U.rows())
// U has the dimension of the source
// since PLU must have dimension equal to the source

//TODO: L matrix could be a lower triangular.

template <typename T, index_t M, index_t N>
class PLUBase {
public:
    SimpleMatrix<T,M,M>  L;
    SimpleMatrix<T,M,N>  U;
    PermutationMatrix<M> P;
    bool singular;
    bool swap_parity;
    
    PLUBase(index_t rows, index_t cols):
        singular(false),
        swap_parity(false) {}
};

template <typename T, index_t M>
class PLUBase<T,M,DYNAMIC_DIM> {
public:
    SimpleMatrix<T,M,M> L;
    SimpleMatrix<T,M,DYNAMIC_DIM> U;
    PermutationMatrix<M> P;
    bool singular;
    bool swap_parity;
    
    PLUBase(index_t rows, index_t cols):
        U(cols),
        singular(false),
        swap_parity(false) {}
};

template <typename T, index_t N>
class PLUBase<T,DYNAMIC_DIM,N> {
public:
    SimpleMatrix<T,DYNAMIC_DIM,DYNAMIC_DIM>  L;
    SimpleMatrix<T,DYNAMIC_DIM,N>  U;
    PermutationMatrix<DYNAMIC_DIM> P;
    bool singular;
    bool swap_parity;
    
    PLUBase(index_t rows, index_t cols):
        L(rows,rows),
        U(rows),
        P(rows),
        singular(false),
        swap_parity(false) {}
};

template <typename T>
class PLUBase<T,DYNAMIC_DIM,DYNAMIC_DIM> {
public:
    SimpleMatrix<T,DYNAMIC_DIM,DYNAMIC_DIM>  L;
    SimpleMatrix<T,DYNAMIC_DIM,DYNAMIC_DIM>  U;
    PermutationMatrix<DYNAMIC_DIM> P;
    bool singular;
    bool swap_parity;
    
    PLUBase(index_t rows, index_t cols):
        L(rows, rows),
        U(rows, cols),
        P(rows),
        singular(false),
        swap_parity(false) {}
};

}; // namespace detail

//////////// PLU class ////////////

template <typename T, index_t M, index_t N>
class PLUDecomposition : public detail::PLUBase<T,M,N> {
    
    template <typename Mx> 
    friend PLUDecomposition<typename Mx::elem_t,Mx::ROWDIM,Mx::COLDIM> plu_decompose(const Mx &src);
    
protected:
    PLUDecomposition(index_t n_r, index_t n_c):
        detail::PLUBase<T,M,N>(n_r, n_c) {}
public:
    
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
            index_t n = this->U.rows();
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
        if ((M == DYNAMIC_DIM or M != K) && this->U.rows() != K) {
            throw DimensionMismatchException(this->U.rows(), 1, K, 1);
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
        if ((J*K == 0 or J != M or K != N) and (into->rows() != this->U.rows() or into->cols() != this->U.cols())) {
            throw DimensionMismatchException(into->rows(), into->cols(), this->U.rows(), this->U.cols());
        }
        #endif
        _inverseTranspose(into->begin());
        into->transpose();
    }
    
    
    T det() const {
        const index_t n = this->L.rows();
        T k = getParity();
        for (index_t i = 0; i < n; i++) {
            k *= this->U.get(i,i);
        }
        return k;
    }
    
    
    inline bool isSingular() const {
        return this->singular;
    }
    
    
    inline int getParity() const {
        return this->swap_parity ? -1 : 1;
    }
    
protected:
    
    inline void _checkIsSquare() const {
        if ((M * N == 0 or M != N) and this->U.rows() != this->U.cols()) {
            throw NonsquareMatrixException(this->U.rows(), this->U.cols());
        }
    }
    
    
    template <typename S>
    void _linearSolve(S *dest, const S *b, bool permute=true) const {
        const index_t  n = this->U.rows();
        const index_t *p = this->P.getRowSources(); // permutation map
        
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
                y[r] -= y[c] * this->L.get(r,c);
            }
        }
        // now with y, we may obtain x from:
        // Ux = y
        for (index_t r = n - 1; r >= 0; r--) {
            // dest[r] = y[r]; // nop; dest and y are the same!
            for (index_t c = n - 1; c > r; c--) {
                dest[r] -= dest[c] * this->U.get(r,c);
            }
            dest[r] /= this->U.get(r,r);
        }
    }
    
    
    template <typename S>
    void _inverseTranspose(S *dest) const {
        // Set LUx = PI and solve for x, choosing columns of PI one at a time.
        // Here, we use the rows of our destination matrix as though they are
        // column vectors, and the caller will transpose.
        const index_t *p = this->P.getColSources();
        const index_t n = this->U.rows();
        std::fill(dest, dest + (n*n), 0);
        for (index_t i = 0; i < n; i++, dest += n) {
            dest[p[i]] = 1;
            _linearSolve(dest, dest, false);
        }
    }
};

//////////// user-facing PLU decompose function ////////////

template <typename Mx>
PLUDecomposition<typename Mx::elem_t, Mx::ROWDIM, Mx::COLDIM> plu_decompose(const Mx& src) {
    
    PLUDecomposition<typename Mx::elem_t, Mx::ROWDIM, Mx::COLDIM> decomp(src.rows(), src.cols());
    
    mtxcopy(&decomp.U, src);
    detail::TemporaryStorage<index_t, Mx::ROWDIM> reorder(src.rows());
    bool ok = detail::_ImplDecompPLU(&decomp.L, &decomp.U, reorder.get(), &decomp.swap_parity);
    
    if (not ok) {
        decomp.singular = true;
    } else {
        decomp.P.setRowSources(reorder.get());
    }
    
    return decomp;
}

//todo:
//plu_decompose(PLUDecomposition*, Mx)
//plu_decompose(*L,*U,*P,Mx)
 // note that must setIdentity() first.

};


#endif /* LUDECOMP_H_ */
