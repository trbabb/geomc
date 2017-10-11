/*
 * PermutationMatrix.h
 *
 *  Created on: Apr 28, 2013
 *      Author: tbabb
 */

#ifndef PERMUTATIONMATRIX_H_
#define PERMUTATIONMATRIX_H_

#include <utility>
#include <geomc/Storage.h>
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
    UniqueStorage<index_t,DYNAMIC_DIM> row_src;
    UniqueStorage<index_t,DYNAMIC_DIM> row_dst;
public:
    explicit PermuteMatrixBase(index_t n):
                row_src(n),
                row_dst(n) {}
    
protected:
    inline index_t _rows() const { return row_src.size(); }
    inline index_t _cols() const { return row_src.size(); }
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
 * Permutation matrices are always (n x n) and have only elements that are zero 
 * or 1. Each row and column has exactly one `1` element.
 * 
 * An (*n x n*) permutation matrix uses _O(n)_ storage, and multiplications by 
 * other matrices are optimized to perform the permutation directly in 
 * _O(n<sup>2</sup>)_ (rather than _O(n<sup>3</sup>)_) time. Multiplication of 
 * two permutation matrices is _O(n)_.
 */
template <index_t N>
class PermutationMatrix : public detail::MatrixBase<bool,N,N, PermutationMatrix<N> >, public detail::PermuteMatrixBase<N> {
private:
    typedef detail::PermuteMatrixBase<N> parent_t;
    
    index_t sign; // 0 means "not computed yet".
                  // could possibly compute on construction.
                  // (or don't store at all).
                  // disadvantage: duplicated code; computeSign() is non-const. :(
    
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
    
    /**
     * Compute the signature of this permutation.
     *
     * Runs in _O(n)_ time; or _O(1)_ time if previously computed.
     *
     * @return -1 if the number of transpositions in this permutation is odd, 1 otherwise.
     */
    index_t getSign() {
        if (sign == 0) {
            index_t n = detail::PermuteMatrixBase<N>::_rows();
            index_t *row_src = parent_t::getSrcData();
            index_t *row_dst = parent_t::getDstData();
            
            sign = computeSign(row_src, row_dst, n) ? -1 : 1;
            
            // we just destroyed our dst data. recompute.
            for (index_t i = 0; i < n; i++) {
                row_dst[row_src[i]] = i;
            }
        } else {
            return sign;
        }
    }
    
    
    /**
     * @return Number of rows in the matrix. Always equal to the number of columns.
     */
    inline index_t rows() const {
        return parent_t::_rows();
    }
    
    
    /**
     * @return Number of columns in the matrix. Always equal to the number of rows.
     */
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
     * return an array `a` such that row `D`<sub>`i`</sub>` = M`<sub>`a[i]`</sub>. 
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
     * return an array `a` such that column `D`<sub>`i`</sub>` = M`<sub>`a[i]`</sub>. 
     */
    inline const index_t *getColSources() const {
        return parent_t::getDstData();
    }
    
    /**
     * Array containing a mapping of source columns to destination columns in a 
     * column-permuting operation. 
     * 
     * In other words, given the column-permuting multiplication:
     * 
     *     M * P = D
     * 
     * return an array `a` such that column `D`<sub>`a[i]`</sub>` = M`<sub>`i`</sub>.
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
     * return an array `a` such that row `D`<sub>`a[i]`</sub>` = M`<sub>`i`</sub>.
     * 
     * Because of the property that a permutation matrix's inverse is its 
     * transpose, this function is equivalent to `getColSources()`.
     */
    inline const index_t *getRowDestinations() const {
        return parent_t::getDstData();
    }
    
    /**
     * Set the permutation described by this matrix by passing a mapping of 
     * destination rows to source rows in a row-permuting operation.
     * 
     * In other words, define the row-permuting multiplicaton:
     * 
     *     P * M = D
     * 
     * with an array `a` such that row `D`<sub>`i`</sub>` = M`<sub>`a[i]`</sub>. 
     * 
     * @param p Array of indecies.
     */
    void setRowSources(const index_t *p) {
        index_t n = detail::PermuteMatrixBase<N>::_rows();
        index_t *row_src = parent_t::getSrcData();
        index_t *row_dst = parent_t::getDstData();
        std::copy(p, p + n, row_src);
        for (index_t i = 0; i < n; i++) {
            row_dst[row_src[i]] = i;
        }
        sign = 0;
    }
    
    /**
     * Set the permutation described by this matrix by passing a mapping of 
     * destination columns to source columns in a column-permuting operation. 
     * 
     * In other words, define the column-permuting multiplicaton:
     * 
     *     M * P = D
     * 
     * with an array `a` such that column `D`<sub>`i`</sub>` = M`<sub>`a[i]`</sub>. 
     * 
     * @param p Array of indecies.
     */
    void setColSources(const index_t *p) {
        index_t n = detail::PermuteMatrixBase<N>::_cols();
        index_t *row_src = parent_t::getSrcData();
        index_t *row_dst = parent_t::getDstData();
        std::copy(p, p + n, row_dst);
        for (index_t i = 0; i < n; i++) {
            row_src[row_dst[i]] = i;
        }
        sign = 0;
    }
    
    /**
     * Set the permutation described by this matrix by passing an array mapping 
     * source rows to destination rows in a row-permuting operation. 
     * 
     * In other words, define the row-permuting multiplication:
     * 
     *     P * M = D
     * 
     * with an array `a` such that row `D`<sub>`a[i]`</sub>` = M`<sub>`i`</sub>.
     * 
     * Because of the property that a permutation matrix's inverse is its 
     * transpose, this function is equivalent to `setColSources()`.
     * 
     * @param p Array of indecies.
     */
    inline void setRowDestinations(const index_t *p) {
        setColSources(p);
    }
    
    /**
     * Set the permutation described by this matrix by passing an array mapping 
     * source columns to destination columns in a column-permuting operation. 
     * 
     * In other words, define the column-permuting multiplication:
     * 
     *     M * P = D
     * 
     * with an array `a` such that column `D`<sub>`a[i]`</sub>` = M`<sub>`i`</sub>.
     * 
     * Because of the property that a permutation matrix's inverse is its 
     * transpose, this function is equivalent to `setRowSources()`.
     * 
     * @param p Array of indecies.
     */
    inline void setColDestinations(const index_t *p) {
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
        sign *= -1;
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
        sign *= -1;
    }
    
    /**
     * Get the element at `(row, col)`.
     * @param row Zero-indexed row coordinate
     * @param col Zero-indexed column coordinate
     * @return The element at `(row, col)`; either 0 or 1.
     */
    inline bool get(index_t row, index_t col) const {
        return parent_t::getSrcData()[row] == col;
    }
    
    /**
     * Reset this matrix to the identity permutation. 
     */
    void setIdentity() {
        index_t *i0 = parent_t::getSrcData();
        index_t *i1 = parent_t::getDstData();
        for (index_t i = 0; i < rows(); i++) {
            *i0++ = i;
            *i1++ = i;
        }
        sign = 1;
    }
    
    /**
     * Transpose this matrix in-place in _O(1)_ time. For this type of matrix, also the inverse matrix.
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
    
private:
    
    // compute the signature of this permutation, by counting the length of the cycles.
    // because we can decompose a permutation into a product of cyclical permutations,
    // we can decompose the sign as the product of the sign of those cycles.
    // a cycle has odd sign (-1) if the number of elements in the cycle is even.
    // of course we do not actually have to keep track of the count; only its parity.
    // note that inverting/transposing a permutation does not change its sign,
    // so we can interchange row_src and row_dst at our convenience without changing the answer.
    // whichever one we are given, we use the other as a buffer to mark which elements
    // have been visited; this allows us to avoid a dynamic memory alloc. for that reason,
    // we must call this function after one index map has been filled, but before the other has been computed.
    // runs in O(n) time.
    
    bool computeSign(index_t* p, index_t* buf, index_t n) {
        std::fill(buf, buf+n, 0);
        
        // follow the cycles.
        bool odd = false;
        for (index_t i = 0; i < n; ++i) {
            if  (buf[i] == 1) continue; // not strictly necessary, but saves us an indirection
            else buf[i] = 1;
            for (index_t j = p[i]; buf[j] != 1; j = p[j]) {
                buf[j] = 1; // mark that this element has been visited.
                odd = !odd; // the product of the sign of the cycles has flipped.
            }
        }
        return odd;
    }
};

};

#endif /* PERMUTATIONMATRIX_H_ */
