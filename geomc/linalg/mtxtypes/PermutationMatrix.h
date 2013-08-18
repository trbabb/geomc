/*
 * PermutationMatrix.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef PERMUTATIONMATRIX_H_
#define PERMUTATIONMATRIX_H_

#include <utility>
#include <geomc/linalg/mtxdetail/MatrixBase.h>

namespace geom {

namespace detail {

// fwd decl for friend function:
template <index_t N>
bool mtxequal(const geom::PermutationMatrix<N> &a, const geom::PermutationMatrix<N> &b);

template <index_t N> 
class PermuteMatrixBase {
    // buffers for indirection table
    index_t data_0[N];
    index_t data_1[N];
    bool swap;
public:
    PermuteMatrixBase():swap(false) {}
    
    explicit PermuteMatrixBase(index_t n):swap(false) {}
    
protected:
    
    inline index_t _rows() const { return N; }
    inline index_t _cols() const { return N; }
          inline index_t *getSrcData()       { return swap ? data_0 : data_1; } // we could avoid this test by 
    const inline index_t *getSrcData() const { return swap ? data_0 : data_1; } // storing pointers, and fixing
          inline index_t *getDstData()       { return swap ? data_1 : data_0; } // them up at object copy to
    const inline index_t *getDstData() const { return swap ? data_1 : data_0; } // prevent aliasing.
    inline void           swapPointers() { swap = not swap; }
};


template <>
class PermuteMatrixBase<DYNAMIC_DIM> {
    boost::shared_array<index_t> row_src;
    boost::shared_array<index_t> row_dst;
    index_t n;
public:
    explicit PermuteMatrixBase(index_t n):
                row_src(new index_t[n]),
                row_dst(new index_t[n]),
                n(n) {}
    
protected:
    inline index_t _rows() const { return n; }
    inline index_t _cols() const { return n; }
    const inline index_t *getSrcData() const { return row_src.get(); }
          inline index_t *getSrcData()       { return row_src.get(); }
    const inline index_t *getDstData() const { return row_dst.get(); } 
          inline index_t *getDstData()       { return row_dst.get(); }
    inline void           swapPointers() { std::swap(row_src, row_dst); }
};


}; // namespace detail

// fwd decl
namespace detail {
template <typename Ma, typename Mb, typename Enable> class _ImplMtxMul;
};

/** @ingroup matrix 
 *  @brief A matrix which, by multiplication, permutes the rows or columns of
 * another matrix.
 * 
 * A permutation matrix `P` permutes rows if left-multiplied (`P * M`), and permutes
 * columns if right-multiplied (`M * P`).
 * 
 * Permutation matrices use _O(n)_ storage, and multiplication operators with them
 * are optimized to perform the permutation directly in _O(n<sup>2</sup>)_
 * (rather than _O(n<sup>3</sup>)_) time.
 * 
 * Permutation matrices are always (n x n) and have only elements that are zero 
 * or 1. Each row and column has exactly one `1` element.
 */
template <index_t N>
class PermutationMatrix : public detail::MatrixBase<bool,N,N, PermutationMatrix<N> >, public detail::PermuteMatrixBase<N> {
private:
    typedef detail::PermuteMatrixBase<N> parent_t;
    
public:
    
    template <typename Md, typename Ma, typename Mb>
    friend class detail::_ImplMtxMul;
    
    /**
     * Construct a new identity permutation matrix.
     */
    PermutationMatrix() {
        setIdentity();
    }
    
    /**
     * Construct a new identity permutation matrix of size `n`.
     * (Dynamic only).
     */
    explicit PermutationMatrix(index_t n):parent_t(n) {
        setIdentity();
    }
    
    //////////// methods ////////////
    
    inline index_t rows() const {
        return parent_t::_rows();
    }
    
    inline index_t cols() const {
        return parent_t::_cols();
    }
    
    /**
     * Array containing a mapping of destination rows to source rows in a 
     * row-permuting operation. 
     * 
     * In other words, given the row-permuting multiplicaton:
     * 
     *     P * M = D
     * 
     * return an array `a` such that row `D[i] = M[a[i]]`. 
     */
    inline const index_t *getRowSources() const {
        return parent_t::getSrcData();
    }
    
    /**
     * Array containing a mapping of destination columns to source columns in a 
     * column-permuting operation. 
     * 
     * In other words, given the column-permuting multiplicaton:
     * 
     *     M * P = D
     * 
     * return an array `a` such that column `D[i] = M[a[i]]`. 
     */
    inline const index_t *getColSources() const {
        return parent_t::getDstData();
    }
    
    void setRowSources(index_t *p) {
        index_t n = detail::PermuteMatrixBase<N>::_rows();
        index_t *row_src = parent_t::getSrcData();
        index_t *row_dst = parent_t::getDstData();
        std::copy(p, p + n, row_src);
        for (index_t i = 0; i < n; i++) {
            row_dst[row_src[i]] = i;
        }
    }
    
    void setColSources(index_t *p) {
        index_t n = detail::PermuteMatrixBase<N>::_cols();
        index_t *row_src = parent_t::getSrcData();
        index_t *row_dst = parent_t::getDstData();
        std::copy(p, p + n, row_dst);
        for (index_t i = 0; i < n; i++) {
            row_src[row_dst[i]] = i;
        }
    }
    
    /**
     * Array containing a mapping of source columns to destination columns in a 
     * column-permuting operation. 
     * 
     * In other words, given the column-permuting multiplication:
     * 
     *     M * P = D
     * 
     * return an array `a` such that column `D[a[i]] = M[i]`.
     * 
     * Because of the property that a permutation matrix's inverse is its 
     * transpose, this function is equivalent to `getRowSources()`.
     */
    inline const index_t *getColDestinations() const {
        return parent_t::getSrcData();
    }
    
    /**
     * Array containing a mapping of source rows to destination rows in a row-permuting
     * operation. 
     * 
     * In other words, given the row-permuting multiplication:
     * 
     *     P * M = D
     * 
     * return an array `a` such that row `D[a[i]] = M[i]`.
     * 
     * Because of the property that a permutation matrix's inverse is its 
     * transpose, this function is equivalent to `getColSources()`.
     */
    inline const index_t *getRowDestinations() const {
        return parent_t::getDstData();
    }
    
    inline void setRowDesintations(index_t *p) {
        setColSources(p);
    }
    
    inline void setColDestinations(index_t *p) {
        setRowSources(p);
    }
    
    ///////////////////////////////////////////////////
    
    /**
     * Adjust this matrix such that rows `a` and `b` are swapped in the destination
     * matrix after applying a row permutation. This operation is cumulative on 
     * any previous swaps.
     * 
     * @param a A row index.
     * @param b A row index.
     */
    void swap_rows(index_t a, index_t b) {
        index_t *src = parent_t::getSrcData();
        index_t *dst = parent_t::getDstData();
        // inverse swap
        dst[src[a]] = b;
        dst[src[b]] = a;
        // swap a/b
        std::swap(src[a], src[b]);
    }
    
    /**
     * Adjust this matrix such that columns `a` and `b` are swapped in the destination
     * matrix after applying a column permutation. This operation is cumulative on 
     * any previous swaps.
     * 
     * @param a A column index.
     * @param b A column index.
     */
    void swap_cols(index_t a, index_t b) {
        index_t *src = parent_t::getSrcData();
        index_t *dst = parent_t::getDstData();
        // inverse swap
        src[dst[a]] = b;
        src[dst[b]] = a;
        // swap a/b
        std::swap(dst[a], dst[b]);
    }
    
    inline bool get(index_t row, index_t col) const {
        return parent_t::getSrcData()[row] == col;
    }
    
    void setIdentity() {
        index_t *i0 = parent_t::getSrcData();
        index_t *i1 = parent_t::getDstData();
        for (index_t i = 0; i < rows(); i++) {
            *i0++ = i;
            *i1++ = i;
        }
    }
    
    /**
     * Matrix transpose. For this type of matrix, also the inverse matrix.
     */
    inline void transpose() {
        detail::PermuteMatrixBase<N>::swapPointers();
    }

    inline void getStorageIDs(storage_id_t *buf) const {
        *buf = this;
    }
    
    inline index_t getStorageIDCount() const {
        return 1;
    }
    
    /******** friends ********/
    
    friend bool detail::mtxequal<N>(const PermutationMatrix<N> &a, const PermutationMatrix<N> &b);
};

};

#endif /* PERMUTATIONMATRIX_H_ */
