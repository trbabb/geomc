#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MatrixIterator

#include <iostream>
#include <random>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <geomc/linalg/Matrix.h>

// because macros can't be made to not eat the comma in
// constructs like MACRO(Thing<A,B>).
// macros are very silly in general and this is silly too
#define SINGLE_ARG(...) __VA_ARGS__

using namespace geom;
using namespace std;

typedef std::mt19937_64 rng_t;


template <typename T, index_t M, index_t N, MatrixLayout L>
void randomize_matrix(SimpleMatrix<T,M,N,L>* mx, rng_t* rng) {
    typedef std::normal_distribution<T> d_normal_t;
    d_normal_t gauss = d_normal_t(0., 1); // (ctr, variance)
    
    for (index_t r = 0; r < mx->rows(); ++r) {
        for (index_t c = 0; c < mx->cols(); ++c) {
            (*mx)(r,c) = gauss(*rng);
        }
    }
}


template <index_t M, index_t N>
struct LayoutFixture {
    SimpleMatrix<float,M,N,ROW_MAJOR> m_a;
    SimpleMatrix<float,M,N,COL_MAJOR> m_b;
    
    LayoutFixture() {
        index_t i = 0;
        for (index_t r = 0; r < M; ++r) {
            for (index_t c = 0; c < N; ++c, ++i) {
                m_a(r,c) = i;
                m_b(r,c) = i;
            }
        }
    }
};


template <typename Iterator>
index_t count_steps(Iterator begin, Iterator end) {
    index_t ct = 0;
    for (Iterator i = begin; i != end; ++i, ++ct) {}
    return ct;
}


namespace std {
    
    template <typename T, bool Const>
    std::ostream& operator<<(std::ostream& o, const geom::TransposeIterator<T,Const>& i) {
        o << "<" << i.base << "[" << (i.p - i.base) << "] = " << (*i) << ">";
        return o;
    }
}


BOOST_FIXTURE_TEST_SUITE(matrix_layout_tests, SINGLE_ARG(LayoutFixture<5,7>))


BOOST_AUTO_TEST_CASE(verify_layout_parallel) {
    // the two matrices were filled using set(r,c).
    // verify that they compare equal.
    BOOST_CHECK(m_a == m_b);
}

BOOST_AUTO_TEST_CASE(verify_matrix_copy) {
    // copy matrix A to matrix B and back; verify equality
    
    for (auto i = m_b.begin(); i != m_b.end(); ++i) {
        // destroy original contents
        *i = 0;
    }
    mtxcopy(&m_b, m_a);
    BOOST_CHECK(m_a == m_b);
    
    for (auto i = m_a.begin(); i != m_a.end(); ++i) {
        // destroy original contents
        *i = 0;
    }
    mtxcopy(&m_a, m_b);
    BOOST_CHECK(m_a == m_b);
}

// xxx fails on end-of-array case
BOOST_AUTO_TEST_CASE(verify_col_iterator_offend) {
    float* base = m_a.data_begin();
    for (index_t i = 0; i < m_a.cols(); ++i) {
        auto col_start = m_a.col(i);
        // the end of a column is the beginning of the next one
        auto col_end   = m_a.col(i + 1);
        float* col_startp = &(*col_start);
        float* col_endp   = &(*col_end);
        // the number of items in a column (i.e. between the start and end ptrs)
        // should equal the number of rows
        BOOST_CHECK_EQUAL(col_end - col_start, m_a.rows());
        // the start of a column 
        BOOST_CHECK_EQUAL(col_start + m_a.rows(), col_end);
    }
}

// xxx fails on the end-of-array case
BOOST_AUTO_TEST_CASE(verify_row_iterator_offend) {
    float* base = m_b.data_begin();
    for (index_t i = 0; i < m_b.rows(); ++i) {
        auto row_start = m_b.row(i);
        // the end of a row is the beginning of the next one
        auto row_end   = m_b.row(i + 1);
        float* row_startp = &(*row_start);
        float* row_endp   = &(*row_end);
        // the number of items in a row (i.e. between the start and end ptrs)
        // should equal the number of cols
        BOOST_CHECK_EQUAL(row_end - row_start, m_b.cols());
        // the start of a row 
        BOOST_CHECK_EQUAL(row_start + m_b.cols(), row_end);
    }
}

BOOST_AUTO_TEST_CASE(verify_iterator_count) {
    // m_a = row major
    // m_b = col major
    index_t ma_sz = m_a.rows() * m_a.cols();
    index_t mb_sz = m_b.rows() * m_b.cols();
    // do the begin/end iterators span the size of the matrix?
    BOOST_CHECK_EQUAL(count_steps(m_a.begin(), m_a.end()), ma_sz);
    BOOST_CHECK_EQUAL(count_steps(m_b.begin(), m_b.end()), mb_sz);
    // do the row iterators span the matrix?
    BOOST_CHECK_EQUAL(count_steps(m_a.row(0), m_a.row(m_a.rows())), ma_sz);
    BOOST_CHECK_EQUAL(count_steps(m_b.row(0), m_b.row(m_b.rows())), mb_sz);
    // do the column iterators span the matrix?
    BOOST_CHECK_EQUAL(count_steps(m_a.col(0), m_a.col(m_a.cols())), ma_sz);
    BOOST_CHECK_EQUAL(count_steps(m_b.col(0), m_b.col(m_b.cols())), mb_sz); 
}

BOOST_AUTO_TEST_CASE(verify_matrix_add) {
    auto m_d = m_a + m_b;
    BOOST_CHECK_EQUAL(m_d.rows(), m_a.rows());
    BOOST_CHECK_EQUAL(m_d.cols(), m_a.cols());
    
    auto i_d = m_d.begin();
    auto i_a = m_a.begin();
    auto i_b = m_b.begin();
    for (; i_d != m_d.end(); ++i_d, ++i_a, ++i_b) {
        BOOST_CHECK_EQUAL(*i_d, *i_a + *i_b);
    }
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(matrix_inv)

// todo: iterate this guy over sizes 2, 3, 4, 5, 6
BOOST_AUTO_TEST_CASE(verify_matrix_inv) {
    rng_t rng(17581355241LL);
    SimpleMatrix<double,3,3,ROW_MAJOR> m0r;
    SimpleMatrix<double,3,3,ROW_MAJOR> m1r;
    SimpleMatrix<double,3,3,COL_MAJOR> m0c;
    SimpleMatrix<double,3,3,COL_MAJOR> m1c;
    
    randomize_matrix(&m0r, &rng);
    
    // verify that inverting into row major and col major matrices is the same.
    inv(&m1r, m0r);
    inv(&m1c, m0r);
    BOOST_CHECK(m1r == m1c);
    
    // verify the inverse is correct, i.e. that M * M^-1 = I
    mul(&m0c, m0r, m1r);
    for (index_t r = 0; r < m0r.rows(); ++r) {
        for (index_t c = 0; c < m0r.cols(); ++c) {
            if (r == c) BOOST_CHECK_CLOSE(m0r(r,c), 1, 0.0001);
            else BOOST_CHECK_SMALL(m0r(r,c), 1e-6);
        }
    }
    
    // now test again with a column-major matrix as the source.
    randomize_matrix(&m0c, &rng);
    inv(&m1c, m0c);
    inv(&m1r, m0c);
    BOOST_CHECK(m1r == m1c);
    
    
    // verify correctness
    mul(&m0r, m0c, m1c);
    for (index_t r = 0; r < m0r.rows(); ++r) {
        for (index_t c = 0; c < m0r.cols(); ++c) {
            if (r == c) BOOST_CHECK_CLOSE(m0r(r,c), 1, 0.0001);
            else BOOST_CHECK_SMALL(m0r(r,c), 1e-6);
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()

