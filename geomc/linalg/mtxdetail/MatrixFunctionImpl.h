/*
 * MatrixFunctionImpl.h
 *
 *  Created on: Jan 28, 2013
 *      Author: tbabb
 */

#ifndef MATRIXFUNCTIONIMPL_H_
#define MATRIXFUNCTIONIMPL_H_

#include <algorithm>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>

/*************************************
 * Matrix copy operations            *
 *************************************/

namespace geom {

    template <typename Ms, typename Md>
    inline void copyMatrixRegion(const Ms &src, Md &dst, const MatrixRegion &src_region, const MatrixCoord &dst_begin){
        std::pair<typename Ms::region_iterator, typename Ms::region_iterator> p = src.region(src_region);
        MatrixRegion dst_region = MatrixRegion(dst_begin, dst_begin + src_region.getDimensions());
        std::copy(p[0], p[1], dst.region(dst_region));
    }
    
};

/*************************************************
 * Matrix instantiators                          *
 *                                               *
 * These are necessary so that all matrix        * 
 * operator wrappers can allocate destination    *
 * matrices without knowing their type or        *
 * constructor arguments.                        *
 *************************************************/

namespace geom {
namespace detail {

template <typename M>
struct _ImplMtxInstance {
    static inline M instance(index_t r, index_t c) {
        return M(r,c);
    }
};

template <index_t N>
struct _ImplMtxInstance <geom::PermutationMatrix<N> > {
    static inline geom::PermutationMatrix<N> instance(index_t r, index_t c) {
        assert(r == c);
        return geom::PermutationMatrix<N>(r);
    }
};

template <typename T, index_t N>
struct _ImplMtxInstance <geom::Vec<T,N> > {
    static inline geom::Vec<T,N> instance(index_t r, index_t c) {
        return geom::Vec<T,N>();
    }
};


/************************************
 * Matrix transpose implementations *
 ************************************/

template <typename Md, typename M>
inline void mtx_copy_txpose(Md *into, const M& m) {
    for (index_t c = 0; c < into->cols(); c++) {
        for (int r = 0; r < into->rows(); r++) {
            into->set(r, c, m.get(c, r));
        }
    }
}

template <typename M> 
struct _ImplMtxTxpose {
    
    typedef geom::SimpleMatrix<typename M::elem_t, M::COLDIM, M::ROWDIM> return_t;

    template <typename Md>
    static void transpose(Md *d, const M& m) {
        mtx_copy_txpose(d, m);
    }
};

template <typename T, index_t M, index_t N>
struct _ImplMtxTxpose<geom::DiagMatrix<T,M,N> > {
    
    typedef geom::DiagMatrix<T,N,M> return_t;
    
    template <typename Md>
    inline static void transpose(Md *d, const geom::DiagMatrix<T,M,N> &m) {
        mtx_copy_txpose(d, m);
    }
    
    inline static void transpose(return_t *d, const geom::DiagMatrix<T,M,N> &m) {
        // diagonal unaffected by transpose
        std::copy(m.diagonal_begin(), m.diagonal_end(), d->diagonal_begin());
    }
};

template <index_t N>
struct _ImplMtxTxpose<geom::PermutationMatrix<N> > {
    
    typedef geom::PermutationMatrix<N> return_t;
    
    template <typename Md>
    inline static void transpose(Md *d, const geom::PermutationMatrix<N> &m) {
        mtx_copy_txpose(d, m);
    }
    
    inline static void transpose(geom::PermutationMatrix<DYNAMIC_DIM> *d, const geom::PermutationMatrix<N> &m) {
        d->setRowSources(m.getColSources()); // better than (*d = m).transpose(), which aliases memory.
    }
    
    template <index_t L>
    inline static void transpose(geom::PermutationMatrix<L> *d, const geom::PermutationMatrix<N> &m) {
        *d = m;
        d->transpose();
    }
};

/************************************
 * Matrix equality                  *
 ************************************/

template <typename Ma, typename Mb>
bool mtxequal(const Ma &a, const Mb &b) {
    typename Ma::const_iterator ai = a.begin();
    typename Mb::const_iterator bi = b.begin();
    
    for (; ai != a.end(); ++ai, ++bi) {
        if (*ai != *bi) {
            return false;
        }
    }
    
    return true;
}

template <typename T, index_t M1, index_t N1, index_t M2, index_t N2>
bool mtxequal(const geom::DiagMatrix<T,M1,N1> &a, const geom::DiagMatrix<T,M2,N2> &b) {
    const T *pa = a.diagonal_begin();
    const T *pb = b.diagonal_begin();
    
    for (; pa != a.diagonal_end(); ++pa, ++pb) {
        if (*pa != *pb) {
            return false;
        }
    }
    
    return true;
}

// also has a fwd decl and friend in PermutationMatrix
template <index_t N>
bool mtxequal(const geom::PermutationMatrix<N> &a, const geom::PermutationMatrix<N> &b) {
    const index_t *a_rsrc = a.getSrcData();
    const index_t *b_rsrc = b.getSrcData();
    
    if (a_rsrc == b_rsrc) return true;
    
    for (index_t i = 0; i < a.rows(); i++) {
        if (a_rsrc[i] != b_rsrc[i]) {
            return false;
        }
    }
    
    return true;
}

/************************************
 * Matrix copy                      *
 ************************************/

template <typename Md, typename Mx>
inline void _mtxcopy(Md *into, const Mx &src) {
    std::copy(src.begin(), src.end(), into->begin()); //also works for vectors. neat!
}


template <index_t N>
inline void _mtxcopy(geom::PermutationMatrix<N> *into, const geom::PermutationMatrix<N> &src) {
    into->setRowSources(src.getRowSources());
}


template <typename Td, typename Ts>
void _mtxcopy(geom::SparseMatrix<Td> *into, const geom::SparseMatrix<Ts> &src) {
    into->setZero();
    typedef typename SparseMatrix<Ts>::nonzero_iterator nz_s;
    for (nz_s i_s = src.nonzero_begin(); i_s != src.nonzero_end(); i_s++) {
        MatrixCoord p = i_s->first;
        into->set(p.row, p.col, i_s->second);
    }
}


// also has a fwd decl and friend in SparseMatrix
template <typename T>
inline void _mtxcopy(geom::SparseMatrix<T> *into, const geom::SparseMatrix<T> &src) {
    into->elements = src.elements;
}


template <typename Td, index_t Md, index_t Nd, 
          typename Ts, index_t Ms, index_t Ns>
void _mtxcopy(geom::DiagMatrix<Td, Md, Nd> *into, const geom::DiagMatrix<Ts, Ms, Ns> &src) {
    std::copy(src.diagonal_begin(), src.diagonal_end(), into->diagonal_begin());
}


template <typename Td, typename Ts, index_t Ms, index_t Ns>
void _mtxcopy(geom::SparseMatrix<Td> *into, const geom::DiagMatrix<Ts, Ms, Ns> &src) {
    into->setZero();
    Ts *p = src.diagonal_begin();
    for (index_t i = 0; p != src.diagonal_end(); i++, p++) {
        into->set(i,i,*p);
    }
}


template <typename Td, index_t N>
void _mtxcopy(geom::SparseMatrix<Td> *into, const geom::PermutationMatrix<N> &src) {
    into->setZero();
    index_t *rowdst = src.getColSources(); 
    for (index_t i = 0; i < src.rows(); i++) {
        into->set(i, rowdst[i], 1);
    }
}


// TODO: a permutation matrix
//       x elem_t == bool
//         x specialized mult for mtxes 
//         x and vecs 
//         - and other permute mtxes
//       x exposure for 
//         x mul(), 
//         x xpose(), 
//         x inv() 
//         x which respect static vs. dynamic nature.
//       x specialized inverse (transpose is inverse)
//         x short-circuit copy if src==dest
//         x make sure the short circuit is not subverted by an alloc/copy within
//           the outer, default exposed function.
//         x basically, make explicit global inv()s and transpose()es for permutationmtxes
//           x this may not be necessary, but you should be sure to skip the aliasing
//             dance for permute matrices, which can invert in-place.
//             > technically this is true of 2x2 matrices as well.
//             > actually, the aliasing handling should be pushed down into the invert
//               implementation, since many/most cases will not need it.
//       x specialized transpose (swap row, col indirection tables)
//       X creation function in PermutedMatrix
//         x PermutedMatrix deprecated
//       - specialized determinant (1 or -1)
//       x flatten() in PermutedMatrix
// TODO: in wrapper, allow case where src == upper_dest


}; /* namespace geom::detail */
}; /* namespace geom */

#endif /* MATRIXFUNCTIONIMPL_H_ */