/*
 * MatrixInv.h
 *
 *  Created on: Jul 18, 2013
 *      Author: tbabb
 */

#ifndef MATRIXINV_H_
#define MATRIXINV_H_

#include <algorithm>

#include <geomc/linalg/LinalgTypes.h>
#include <geomc/linalg/LUDecomp.h>

namespace geom {
namespace detail {


#define DET2x2(a,b,c,d) ((a)*(d) - (c)*(b))

/*************************************************
 * General matrix inversion                      *
 *                                               *
 * Inverts the data in an iterator-like object   *
 *************************************************/

//TODO: These are templated over iterator type, which is good
//      in that we don't build redundant inversion functions
//      for all matrices that have raw pointer storage. However,
//      it would be even better if we could separate out the copy
//      step and avoid re-instantiating the entire function for all the
//      myriad matrices that have specialized iterators.

//TODO: for N=3,4 and T in {float,double}, we could implement SSE-optimized 
//      versions of the matrix inverters. Unclear if that would beat -O3.


// NxN static inverse spec
template <typename T, index_t N>
struct _ImplMtxInvRaw {
    
    template <typename S, index_t J, index_t K, StoragePolicy P, typename Mx>
    static inline bool inv(SimpleMatrix<S,J,K,P> *into, const Mx &m) {
        PLUDecomposition<typename Mx::elem_t, Mx::ROWDIM, Mx::COLDIM> plu(m);
        if (plu.isSingular()) return false;
        plu.inverse(into);
        return true;
    }
    
};

// dynamic matrix inverse
template <typename T>
struct _ImplMtxInvRaw<T,DYNAMIC_DIM> {
    
    template <typename Md, typename Mx>
    static bool inv(Md *into, const Mx &m) {
        switch(m.rows()) {
            case 2:
                return geom::detail::_ImplMtxInvRaw<T,2>::inv(into, m);
            case 3:
                return geom::detail::_ImplMtxInvRaw<T,3>::inv(into, m);
            case 4:
                return geom::detail::_ImplMtxInvRaw<T,4>::inv(into, m);
            default:
                PLUDecomposition<typename Mx::elem_t, Mx::ROWDIM, Mx::COLDIM> plu(m);
                if (plu.isSingular()) return false;
                plu.inverse(into);
                return true;
        }
    }
};

// 4x4 inverse spec
template <typename T>
struct _ImplMtxInvRaw<T,4>{
    
    // this allows us to template _inv() for a particular T /only/ over the
    // output iterator, saving some code bloat. also keeps our matrix
    // elements syntactically accessible.
    struct _buf {
        T a00, a01, a02, a03,
          a10, a11, a12, a13,
          a20, a21, a22, a23,
          a30, a31, a32, a33;
    };
    
    template <typename Md, typename Mx>
    static bool inv(Md *into, const Mx &m) {
        _buf b;
        std::copy(m.begin(), m.end(), &b.a00);
        if (_inv(into->begin(), b)) {
            return true;
        }
        return false;
    }
    
    // prevent a little code bloat by templating over iterator type.
    // this allows all the "contiguous memory" matrices to share the
    // same inversion function. Would save even more bloat if we could
    // assign back to <s> inside _inv(), but sadly c++ is dumb
    // and doesn't allow the array initialization syntax for assignments,
    // thus mandating two copies instead of one (one inside the function,
    // and one outside).
    
    template <typename iter_out>
    static bool _inv(iter_out into, const _buf &s) {
        T s0 = DET2x2(s.a00, s.a01, s.a10, s.a11);
        T s1 = DET2x2(s.a00, s.a02, s.a10, s.a12);
        T s2 = DET2x2(s.a00, s.a03, s.a10, s.a13);
        T s3 = DET2x2(s.a01, s.a02, s.a11, s.a12);
        T s4 = DET2x2(s.a01, s.a03, s.a11, s.a13);
        T s5 = DET2x2(s.a02, s.a03, s.a12, s.a13);
        
        T c0 = DET2x2(s.a20, s.a21, s.a30, s.a31);
        T c1 = DET2x2(s.a20, s.a22, s.a30, s.a32);
        T c2 = DET2x2(s.a20, s.a23, s.a30, s.a33);
        T c3 = DET2x2(s.a21, s.a22, s.a31, s.a32);
        T c4 = DET2x2(s.a21, s.a23, s.a31, s.a33);
        T c5 = DET2x2(s.a22, s.a23, s.a32, s.a33);
        
        T det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0;
        
        if (det == 0) return false;
        
        T inv[4][4] = {
                {( s.a11*c5 - s.a12*c4 + s.a13*c3) / det,
                 (-s.a01*c5 + s.a02*c4 - s.a03*c3) / det,
                 ( s.a31*s5 - s.a32*s4 + s.a33*s3) / det,
                 (-s.a21*s5 + s.a22*s4 - s.a23*s3) / det},
                 
                {(-s.a10*c5 + s.a12*c2 - s.a13*c1) / det,
                 ( s.a00*c5 - s.a02*c2 + s.a03*c1) / det,
                 (-s.a30*s5 + s.a32*s2 - s.a33*s1) / det,
                 ( s.a20*s5 - s.a22*s2 + s.a23*s1) / det},
                 
                {( s.a10*c4 - s.a11*c2 + s.a13*c0) / det,
                 (-s.a00*c4 + s.a01*c2 - s.a03*c0) / det,
                 ( s.a30*s4 - s.a31*s2 + s.a33*s0) / det,
                 (-s.a20*s4 + s.a21*s2 - s.a23*s0) / det},
                 
                {(-s.a10*c3 + s.a11*c1 - s.a12*c0) / det,
                 ( s.a00*c3 - s.a01*c1 + s.a02*c0) / det,
                 (-s.a30*s3 + s.a31*s1 - s.a32*s0) / det,
                 ( s.a20*s3 - s.a21*s1 + s.a22*s0) / det}
        };
        
        T *inv_ptr = inv[0];
        std::copy(inv_ptr, inv_ptr+16, into);
        return true;
    }
};


// 3x3 inverse spec
template <typename T>
struct _ImplMtxInvRaw<T,3> {
    
    struct _buf {
        T a, b, c,
          d, e, f,
          g, h, i;
    };
    
    template <typename Md, typename Mx>
    static bool inv(Md *into, const Mx &m) {
        _buf b;
        std::copy(m.begin(), m.end(), &b.a);
        return _inv(into->begin(), b);
    }
    
    template <typename iter_out>
    static bool _inv(iter_out into, const _buf &s) {
        // inverse 3x3 matrix by cramer's rule
        T det = s.a*(s.i*s.e - s.h*s.f) - s.d*(s.i*s.b - s.h*s.c) + s.g*(s.f*s.b - s.e*s.c);
        if (det == 0) return false;
        T inv[3][3] = {{ (s.i*s.e - s.h*s.f)/det, (s.h*s.c - s.i*s.b)/det, (s.f*s.b - s.e*s.c)/det },
                       { (s.g*s.f - s.i*s.d)/det, (s.i*s.a - s.g*s.c)/det, (s.d*s.c - s.f*s.a)/det },
                       { (s.h*s.d - s.g*s.e)/det, (s.g*s.b - s.h*s.a)/det, (s.e*s.a - s.d*s.b)/det }};
        T *from = &(inv[0][0]);
        std::copy(from, from+9, into);
        return true;
    }
};


// 2x2 inverse spec
template <typename T>
struct _ImplMtxInvRaw<T,2>{
    
    template <typename Md, typename Mx>
    static bool inv(Md *into, const Mx &m) {
        // inverse 2x2 matrix by cramer's rule
        T a = m.get(0,0);
        T b = m.get(0,1);
        T c = m.get(1,0);
        T d = m.get(1,1);
        return _inv(into->begin(), a, b, c, d);
    }
    
    template <typename iter_out, typename iter_in>
    static inline bool _inv(iter_out into, iter_in begin, iter_in end) {
        T a[4];
        std::copy(begin, end, a);
        return _inv(into, a[0], a[1], a[2], a[3]);
    }
    
    template <typename iter_out>
    static bool _inv(iter_out into, T a, T b, T c, T d) {
        T det = a*d - c*b;
        if (det == 0) return false;
        T inv[2][2] = {{ d/det, -b/det},
                       {-c/det,  a/det}};
        T *from = inv[0];
        std::copy(from, from+4, into);
        return true;
    }
};


/*************************************
 * Matrix class inversion            *
 *************************************/

//TODO: for icc compiler with floats, use:
//        ftp://download.intel.com/design/PentiumIII/sml/24504301.pdf
//      or:
//        https://github.com/LiraNuna/glsl-sse2/blob/master/source/mat4.h#L324

template <typename Mx>
class _ImplMtxInv {
public:
    
    typedef geom::SimpleMatrix<typename Mx::elem_t, 
                              (Mx::ROWDIM != DYNAMIC_DIM ? Mx::ROWDIM : Mx::COLDIM), 
                              (Mx::ROWDIM != DYNAMIC_DIM ? Mx::ROWDIM : Mx::COLDIM)> return_t;
    
    template <typename Md>
    static inline bool inv(Md *into, const Mx& m) {
        return geom::detail::_ImplMtxInvRaw<typename Md::elem_t, (Mx::ROWDIM != DYNAMIC_DIM ? Mx::ROWDIM : Mx::COLDIM)>::inv(into, m);
    }
};

template <typename T, index_t M, index_t N>
class _ImplMtxInv<geom::DiagMatrix<T,M,N> > {
public:
    
    typedef geom::DiagMatrix<T, (M != 0 ? M : N), (M != 0 ? M : N)> return_t;
    
    template <typename Md>
    static bool inv(Md *into, const geom::DiagMatrix<T,M,N>& m) {
        into->setZero(); // for newly allocated matrices, redundant. oh well.
        T *p = m.diagonal_begin();
        for (index_t i = 0; p != m.diagonal_end(); ++i, ++p) {
            T elem = *p;
            if (elem == 0) return false;
            into->set(i,i, 1/elem);
        }
        return true;
    }
};

// inverse of permutation matrix is its transpose.
template <index_t N>
class _ImplMtxInv<geom::PermutationMatrix<N> > {
    typedef geom::PermutationMatrix<N> M;
public:
    
    typedef geom::PermutationMatrix<N> return_t;
    
    template <typename Md>
    static inline bool inv(Md *into, const M& m) {
        _ImplMtxTxpose<M>::transpose(into, m);
        return true;
    }
};

// workaround likely clang bug.
// clang tries to make a bogus template substitution (Matrix for bool)
// and chokes.
template <typename T>
class _ImplMtxInv<T*> {
    
};

#undef DET2x2

}; // namespace detail
}; // namespace geom

#endif /* MATRIXINV_H_ */
