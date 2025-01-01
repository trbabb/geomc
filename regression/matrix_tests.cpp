#define TEST_MODULE_NAME MatrixTests

#include <iostream>
#include <random>
#include <gtest/gtest.h>
#include <geomc/linalg/Matrix.h>

// xxx: todo: unfuck the matrix inversion code
// xxx: todo: make begin() and end() iterate over memory, not rows.
//            if you want row or column iterators specifically, use row() or col().

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


template <typename MxA, typename MxB>
void check_matrices_close(const MxA& mxa, const MxB& mxb) {
    EXPECT_EQ(mxa.rows(), mxb.rows());
    EXPECT_EQ(mxa.cols(), mxb.cols());
    for (index_t r = 0; r < mxa.rows(); ++r) {
        for (index_t c = 0; c < mxa.cols(); ++c) {
            EXPECT_NEAR(mxa(r,c), mxb(r,c), 0.0001);
        }
    }
}


template <index_t M, index_t N>
struct LayoutFixture : public ::testing::Test {
    typedef SimpleMatrix<float,M,N,ROW_MAJOR> rmaj_t;
    typedef SimpleMatrix<float,M,N,COL_MAJOR> cmaj_t;
    rmaj_t m_a;
    cmaj_t m_b;
    
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

using LayoutFixture57 = LayoutFixture<5,7>;


TEST_F(LayoutFixture57, verify_layout_parallel) {
    // the two matrices were filled using set(r,c).
    // verify that they compare equal.
    EXPECT_TRUE(m_a == m_b);
}

TEST_F(LayoutFixture57, verify_matrix_copy) {
    // copy matrix A to matrix B and back; verify equality
    
    for (auto i = m_b.begin(); i != m_b.end(); ++i) {
        // destroy original contents
        *i = 0;
    }
    mtxcopy(&m_b, m_a);
    EXPECT_TRUE(m_a == m_b);
    
    for (auto i = m_a.begin(); i != m_a.end(); ++i) {
        // destroy original contents
        *i = 0;
    }
    mtxcopy(&m_a, m_b);
    EXPECT_TRUE(m_a == m_b);
}

TEST_F(LayoutFixture57, verify_matrix_copy_init) {
    float x = 3.2;
    for (auto i = m_a.begin(); i != m_a.end(); ++i, x *= 1.07) {
        *i = x;
    }
    cmaj_t m_c(m_a);
    EXPECT_TRUE(m_a == m_c);
}

// xxx fails on end-of-array case
TEST_F(LayoutFixture57, verify_col_iterator_offend) {
    for (index_t i = 0; i < m_a.cols(); ++i) {
        auto col_start = m_a.col(i);
        // the end of a column is the beginning of the next one
        auto col_end   = m_a.col(i + 1);
        // the number of items in a column (i.e. between the start and end ptrs)
        // should equal the number of rows
        EXPECT_EQ(col_end - col_start, m_a.rows());
        // the start of a column 
        EXPECT_EQ(col_start + m_a.rows(), col_end);
    }
}

TEST_F(LayoutFixture57, verify_row_iterator_offend) {
    for (index_t i = 0; i < m_b.rows(); ++i) {
        auto row_start = m_b.row(i);
        // the end of a row is the beginning of the next one
        auto row_end   = m_b.row(i + 1);
        // the number of items in a row (i.e. between the start and end ptrs)
        // should equal the number of cols
        EXPECT_EQ(row_end - row_start, m_b.cols());
        // the start of a row 
        EXPECT_EQ(row_start + m_b.cols(), row_end);
    }
}

TEST_F(LayoutFixture57, verify_iterator_count) {
    // m_a = row major
    // m_b = col major
    index_t ma_sz = m_a.rows() * m_a.cols();
    index_t mb_sz = m_b.rows() * m_b.cols();
    // do the begin/end iterators span the size of the matrix?
    EXPECT_EQ(count_steps(m_a.begin(), m_a.end()), ma_sz);
    EXPECT_EQ(count_steps(m_b.begin(), m_b.end()), mb_sz);
    // do the row iterators span the matrix?
    EXPECT_EQ(count_steps(m_a.row(0), m_a.row(m_a.rows())), ma_sz);
    EXPECT_EQ(count_steps(m_b.row(0), m_b.row(m_b.rows())), mb_sz);
    // do the column iterators span the matrix?
    EXPECT_EQ(count_steps(m_a.col(0), m_a.col(m_a.cols())), ma_sz);
    EXPECT_EQ(count_steps(m_b.col(0), m_b.col(m_b.cols())), mb_sz); 
}

TEST_F(LayoutFixture57, verify_matrix_add) {
    auto m_d = m_a + m_b;
    EXPECT_EQ(m_d.rows(), m_a.rows());
    EXPECT_EQ(m_d.cols(), m_a.cols());
    
    auto i_d = m_d.begin();
    auto i_a = m_a.begin();
    auto i_b = m_b.begin();
    for (; i_d != m_d.end(); ++i_d, ++i_a, ++i_b) {
        EXPECT_EQ(*i_d, *i_a + *i_b);
    }
}

// BOOST_AUTO_TEST_CASE(matrix_speed) {
//     rng_t rng(17581355241LL);
//     static constexpr index_t K = 512;
//     static constexpr index_t L = 1 << 14;
//     using Mx = SimpleMatrix<double,24,24>;
//     Mx* mx = new Mx[K];
//     for (index_t i = 0 ; i < K; ++i) {
//         randomize_matrix(mx + i, &rng);
//     }
    
//     Mx into;
//     double m = 0;
//     auto t0 = std::chrono::high_resolution_clock::now();
//     for (index_t i = 0; i < L; ++i) {
//         inv(&into, mx[i % K]);
//         m = fmod(m + into(0,0), 1);
//     }
//     auto t1 = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> dur = t1 - t0;
//     auto secs = dur.count();
//     std::cout << into << "\n";
//     std::cout << "ops/sec: " << L / secs << "\n";
// }

template <typename T, index_t N>
void test_inv(rng_t& rng) {
    SimpleMatrix<T,N,N,ROW_MAJOR> m0r;
    SimpleMatrix<T,N,N,ROW_MAJOR> m1r;
    SimpleMatrix<T,N,N,COL_MAJOR> m0c;
    SimpleMatrix<T,N,N,COL_MAJOR> m1c;
    SimpleMatrix<T,N,N,ROW_MAJOR> m1d;
    SimpleMatrix<T,N,N,ROW_MAJOR> m1d_buf;
    
    randomize_matrix(&m0r, &rng);
    mtxcopy(&m1d_buf, m0r);
    
    // verify that the inversion succeeds
    EXPECT_TRUE(inv(&m1r, m0r));
    EXPECT_TRUE(inv(&m1c, m0r));
    // check "static" inversion against dynamic inversion:
    EXPECT_TRUE(invNxN(m1d.data_begin(), m1d_buf.data_begin(), N));
    
    // verify that all the methods of inverting produce (nearly) identical results
    check_matrices_close(m1r, m1c);
    check_matrices_close(m1d, m1c);
    
    // verify the inverse is correct, i.e. that M * M^-1 = I
    mul(&m0c, m0r, m1r);
    for (index_t r = 0; r < m0r.rows(); ++r) {
        for (index_t c = 0; c < m0r.cols(); ++c) {
            if (r == c) EXPECT_NEAR(m0c(r,c), 1, 0.0001);
            else EXPECT_NEAR(m0c(r,c), 0, 1e-6);
        }
    }
    
    // now test again with a column-major matrix as the source.
    randomize_matrix(&m0c, &rng);
    inv(&m1c, m0c);
    inv(&m1r, m0c);
    check_matrices_close(m1r, m1c);
    
    // verify correctness
    mul(&m0r, m0c, m1c);
    for (index_t r = 0; r < m0r.rows(); ++r) {
        for (index_t c = 0; c < m0r.cols(); ++c) {
            if (r == c) EXPECT_NEAR(m0r(r,c), 1, 0.0001);
            else EXPECT_NEAR(m0r(r,c), 0, 1e-6);
        }
    }
}

TEST(TEST_MODULE_NAME, verify_matrix_inv) {
    rng_t rng(17581355241LL);
    test_inv<double,2>(rng);
    test_inv<double,3>(rng);
    test_inv<double,4>(rng);
    test_inv<double,5>(rng);
    test_inv<double,6>(rng);
}
