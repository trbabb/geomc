#pragma once
/*
 * MatrixFunctionImpl.h
 *
 *  Created on: Jan 28, 2013
 *      Author: tbabb
 */

#include <algorithm>
#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/mtxdetail/MatrixGlue.h>

/*************************************
 * Matrix copy operations            *
 *************************************/

namespace geom {

template <typename Ms, typename Md>
inline void copyMatrixRegion(const Ms& src, Md& dst, const MatrixRegion& src_region, const Vec<index_t,2>& dst_begin){
    std::pair<typename Ms::region_iterator, typename Ms::region_iterator> p = src.region(src_region);
    MatrixRegion dst_region = MatrixRegion(dst_begin, dst_begin + src_region.dimensions());
    std::copy(p[0], p[1], dst.region(dst_region));
}

/// Transpose a square matrix
template <typename T>
void transpose_square_matrix(T* m, index_t n) {
    detail::MxWrap<T,true> mx = {m, n, n};
    for (index_t i = 0; i < n; ++i) {
        for (index_t j = i + 1; j < n; ++j) {
            std::swap(mx(i,j), mx(j,i));
        }
    }
}


/// Transpose a square matrix with static dimensions
template <typename T, index_t N>
void transpose_square_matrix(T* m) {
    detail::MxWrap<T,true> mx = {m, N, N};
    for (index_t i = 0; i < N; ++i) {
        for (index_t j = i + 1; j < N; ++j) {
            std::swap(mx(i,j), mx(j,i));
        }
    }
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
inline void mtx_copy_txpose(Md* into, const M& m) {
    for (index_t c = 0; c < into->cols(); c++) {
        for (index_t r = 0; r < into->rows(); r++) {
            into->set(r, c, m(c, r));
        }
    }
}

template <typename M> 
struct _ImplMtxTxpose {
    
    typedef geom::SimpleMatrix<typename M::elem_t, M::COLDIM, M::ROWDIM> return_t;

    template <typename Md>
    static void transpose(Md* d, const M& m) {
        mtx_copy_txpose(d, m);
    }
};

template <typename T, index_t M, index_t N, MatrixLayout Lyt, StoragePolicy P>
struct _ImplMtxTxpose<geom::SimpleMatrix<T,M,N,Lyt,P> > {
    
    typedef geom::SimpleMatrix<T,M,N,Lyt,P> Mx;
    typedef geom::SimpleMatrix<T,N,M,Lyt,P> return_t;
    
    
    // xxx: this is masked by alias checking, which will duplicate
    //      the src matrix before we even ask this implementation
    //      how to transpose. this conspires to ensure `d != &m`.
    //      this is instead handled in a specialized function in Matrix.h.
    // In-place transpose a matrix, if we can.
    // This is possible iff the matrix is square.
    // static void transpose(return_t* d, Mx& m) {
    //     if (&m == d and m.rows() == m.cols()) {
    //         for (index_t r = 0; r < m.rows(); r++) {
    //             for (index_t c = 0; c < r; c++) {
    //                 std::swap(m(r,c), m(c,r));
    //             }
    //         }
    //     } else {
    //         // not in-place; fall back
    //         mtx_copy_txpose(d, m);
    //     }
    // }
    
    // not the same type of matrix; therefore not the same matrix.
    // fall back to ordinary copy-transpose.
    template <typename Md>
    inline static void transpose(Md* d, const Mx& m) {
        mtx_copy_txpose(d, m);
    }
};

template <typename T, index_t M, index_t N>
struct _ImplMtxTxpose<geom::DiagMatrix<T,M,N> > {
    
    typedef geom::DiagMatrix<T,N,M> return_t;
    
    template <typename Md>
    inline static void transpose(Md* d, const geom::DiagMatrix<T,M,N>& m) {
        mtx_copy_txpose(d, m);
    }
    
    inline static void transpose(return_t* d, const geom::DiagMatrix<T,M,N>& m) {
        // diagonal unaffected by transpose
        std::copy(m.diagonal_begin(), m.diagonal_end(), d->diagonal_begin());
    }
};

template <index_t N>
struct _ImplMtxTxpose<geom::PermutationMatrix<N> > {
    
    typedef geom::PermutationMatrix<N> return_t;
    
    template <typename Md>
    inline static void transpose(Md* d, const geom::PermutationMatrix<N>& m) {
        mtx_copy_txpose(d, m);
    }
    
    inline static void transpose(geom::PermutationMatrix<DYNAMIC_DIM>* d, const geom::PermutationMatrix<N>& m) {
        d->setRowSources(m.getColSources()); // better than (*d = m).transpose(), which aliases memory.
    }
    
    template <index_t L>
    inline static void transpose(geom::PermutationMatrix<L>* d, const geom::PermutationMatrix<N>& m) {
        *d = m;
        d->transpose();
    }
};

/************************************
 * Matrix equality                  *
 ************************************/

template <typename Ma, typename Mb>
bool mtxequal(const Ma& a, const Mb& b) {
    typename Ma::const_iterator ai = a.begin();
    typename Mb::const_iterator bi = b.begin();
    
    for (; ai != a.end(); ++ai, ++bi) {
        if (*ai != *bi) {
            return false;
        }
    }
    
    return true;
}

// two flat matrices with the same layout can compare bare arrays
template <typename T, typename S, 
          index_t M0, index_t N0, 
          index_t M1, index_t N1, 
          MatrixLayout Lyt, 
          StoragePolicy P0, StoragePolicy P1>
bool mtxequal(
        const SimpleMatrix<T,M0,N0,Lyt,P0>& a,
        const SimpleMatrix<T,M1,N1,Lyt,P1>& b) {
    const T* ai = a.data_begin();
    const T* bi = b.data_begin();
    for (; ai != a.data_end(); ++ai, ++bi) {
        if (*ai != *bi) return false;
    }
    return true;
}

template <typename T, index_t M1, index_t N1, index_t M2, index_t N2>
bool mtxequal(const geom::DiagMatrix<T,M1,N1>& a, const geom::DiagMatrix<T,M2,N2>& b) {
    const T* pa = a.diagonal_begin();
    const T* pb = b.diagonal_begin();
    
    for (; pa != a.diagonal_end(); ++pa, ++pb) {
        if (*pa != *pb) {
            return false;
        }
    }
    
    return true;
}

// also has a fwd decl and friend in PermutationMatrix
template <index_t N>
bool mtxequal(const geom::PermutationMatrix<N>& a, const geom::PermutationMatrix<N>& b) {
    const index_t* a_rsrc = a.getSrcData();
    const index_t* b_rsrc = b.getSrcData();
    
    if (a_rsrc == b_rsrc) return true;
    
    for (index_t i = 0; i < a.rows(); i++) {
        if (a_rsrc[i] != b_rsrc[i]) {
            return false;
        }
    }
    
    return true;
}


// TODO: in wrapper, allow case where src == upper_dest


}; /* namespace geom::detail */
}; /* namespace geom */
