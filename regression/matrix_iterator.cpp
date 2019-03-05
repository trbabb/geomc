#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MatrixIterator

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <geomc/linalg/mtxdetail/MatrixLayout.h>

using namespace geom;
using namespace std;

namespace std {
    template <typename T, bool Const>
    std::ostream& operator<<(std::ostream& o, const geom::StrideIterator<T,Const>& i) {
        o << "<" << i.base << "[" << i.i << "]>";
        return o;
    }
}

#define TEST_ROWS 15
#define TEST_COLS  7


// xxx todo: test negative offsets


typedef FlatMatrixLayout<index_t, ROW_MAJOR> rowmaj;
typedef FlatMatrixLayout<index_t, COL_MAJOR> colmaj;


struct MatrixFixture {
    const index_t rows;
    const index_t cols;
    index_t rowmtx[TEST_ROWS * TEST_COLS];
    index_t colmtx[TEST_ROWS * TEST_COLS];
    
    MatrixFixture():
            rows(TEST_ROWS),
            cols(TEST_COLS) {
        index_t i = 0;
        for (index_t r = 0; r < rows; ++r) {
            for (index_t c = 0; c < cols; ++c, ++i) {
                rowmtx[r * cols + c] = i;
                colmtx[c * rows + r] = i;
            }
        }
    }
    
    virtual ~MatrixFixture() {}
};


BOOST_FIXTURE_TEST_SUITE(matrix_iterator, MatrixFixture)


BOOST_AUTO_TEST_CASE(verify_mtx_index) {
    // verify that row and column layouts are mutual transposes
    for (index_t r = 0; r < rows; ++r) {
        for (index_t c = 0; c < cols; ++c) {
            BOOST_CHECK_EQUAL(
                rowmaj::index(r,c, rows, cols), 
                colmaj::index(c,r, cols, rows));
        }
    }
}


BOOST_AUTO_TEST_CASE(verify_row_order) {
    // verify that row and column layouts agree about the ordering of/within individual rows
    index_t i = 0;
    for (index_t r = 0; r < rows; ++r) {
        auto r_i = rowmaj::row<false>(rowmtx, r, rows, cols);
        auto c_i = colmaj::row<false>(colmtx, r, rows, cols);
        for (index_t c = 0; c < cols; ++c, ++r_i, ++c_i, ++i) {
            BOOST_CHECK_EQUAL(*r_i, *c_i); // agreement?
            BOOST_CHECK_EQUAL(*r_i, i);    // not just wrong in the same way?
        }
    }
}

BOOST_AUTO_TEST_CASE(verify_col_order) {
    // verify that row and column layouts agree about the ordering of/within individual columns
    for (index_t c = 0; c < cols; ++c) {
        auto r_i = rowmaj::col<false>(rowmtx, c, rows, cols);
        auto c_i = colmaj::col<false>(colmtx, c, rows, cols);
        for (index_t r = 0; r < rows; ++r, ++r_i, ++c_i) {
            BOOST_CHECK_EQUAL(*r_i, *c_i);
        }
    }
}

BOOST_AUTO_TEST_CASE(verify_row_iteration) {
    // verify that row and column layouts row-traverse the entire array in the same order
    // (may fail, e.g., at row or column boundaries)
    auto r_ir = rowmaj::row<false>(rowmtx, 0, rows, cols);
    auto c_ir = colmaj::row<false>(colmtx, 0, rows, cols);
    for (index_t i = 0; i < rows * cols; ++i, ++r_ir, ++c_ir) {
        BOOST_CHECK_EQUAL(*r_ir, *c_ir);
    }
}


BOOST_AUTO_TEST_CASE(verify_col_iteration) {
    // verify that row and column layouts column-traverse the entire array in the same order
    // (may fail, e.g., at row or column boundaries)
    auto r_ic = rowmaj::col<false>(rowmtx, 0, rows, cols);
    auto c_ic = colmaj::col<false>(colmtx, 0, rows, cols);
    for (index_t i = 0; i < rows * cols; ++i, ++r_ic, ++c_ic) {
        BOOST_CHECK_EQUAL(*r_ic, *c_ic);
    }
}


BOOST_AUTO_TEST_CASE(verify_offset_increment_row) {
    // verify that "offsetting by `n`" and "incrementing `n` times" are the same for row iterators.
    for (index_t r = 0; r < rows; ++r) {
        auto r_i = rowmaj::row<false>(rowmtx, r, rows, cols);
        auto c_i = colmaj::row<false>(colmtx, r, rows, cols);
        for (index_t c = 0; c < cols; ++c, ++r_i, ++c_i) {
            auto r_j = rowmaj::row<false>(rowmtx, r, rows, cols) + c;
            auto c_j = colmaj::row<false>(colmtx, r, rows, cols) + c;
            // iters consider themselves equal?
            BOOST_CHECK_EQUAL(r_i, r_j);
            BOOST_CHECK_EQUAL(c_i, c_j);
            // iters point to same location?
            BOOST_CHECK_EQUAL(*r_i, *r_j);
            BOOST_CHECK_EQUAL(*c_i, *c_j);
        }
    }
}


BOOST_AUTO_TEST_CASE(verify_offset_increment_col) {
    // verify that "offsetting by `n`" and "incrementing `n` times" are the same for column iterators.
    for (index_t c = 0; c < cols; ++c) {
        auto r_i = rowmaj::col<false>(rowmtx, c, rows, cols);
        auto c_i = colmaj::col<false>(colmtx, c, rows, cols);
        for (index_t r = 0; r < rows; ++r, ++r_i, ++c_i) {
            auto r_j = rowmaj::col<false>(rowmtx, c, rows, cols) + r;
            auto c_j = colmaj::col<false>(colmtx, c, rows, cols) + r;
            // iters consider themselves equal?
            BOOST_CHECK_EQUAL(r_i, r_j);
            BOOST_CHECK_EQUAL(c_i, c_j);
            // iters point to same location?
            BOOST_CHECK_EQUAL(*r_i, *r_j);
            BOOST_CHECK_EQUAL(*c_i, *c_j);
        }
    }
}


BOOST_AUTO_TEST_CASE(verify_distance) {
    // verify that distance() and advance() agree
    for (index_t r_src = 0; r_src < rows; ++r_src) {
        for (index_t c_src = 0; c_src < cols; ++c_src) {
            auto src_i_r = rowmaj::row<false>(rowmtx, r_src, rows, cols) + c_src;
            auto src_i_c = colmaj::row<false>(colmtx, r_src, rows, cols) + c_src;
            for (index_t r_dst = 0; r_dst < rows; ++r_dst) {
                for (index_t c_dst = 0; c_dst < cols; ++c_dst) {
                    auto dst_i_r = rowmaj::row<false>(rowmtx, r_dst, rows, cols) + c_dst;
                    auto dst_i_c = colmaj::row<false>(colmtx, r_dst, rows, cols) + c_dst;
                    auto r_dx = dst_i_r - src_i_r;
                    auto c_dx = dst_i_c - src_i_c;
                    // iterators consider themselves the same?
                    BOOST_CHECK_EQUAL(src_i_r + r_dx, dst_i_r);
                    BOOST_CHECK_EQUAL(src_i_c + c_dx, dst_i_c);
                    // iters point to same location?
                    BOOST_CHECK_EQUAL(*(src_i_r + r_dx), *(dst_i_r));
                    BOOST_CHECK_EQUAL(*(src_i_c + c_dx), *(dst_i_c));
                }
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()
