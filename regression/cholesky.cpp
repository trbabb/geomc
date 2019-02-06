#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Cholesky

// #include <iostream>

#include <random>
#include <boost/test/unit_test.hpp>
#include <geomc/linalg/Cholesky.h>

using namespace geom;
using namespace std;

typedef std::mt19937_64 rng_t;


template <typename T>
SimpleMatrix<T,0,0> random_matrix(index_t sz, rng_t* rng) {
    typedef std::normal_distribution<T> d_normal_t;
    d_normal_t N = d_normal_t(0., 1); // (ctr, variance)
    SimpleMatrix<T,0,0> mx(sz,sz);
    
    for (index_t r = 0; r < sz; ++r) {
        for (index_t c = 0; c < sz; ++c) {
            mx[r][c] = N(*rng);
        }
    }
    return mx;
}


template <typename T>
T matrix_diff(const SimpleMatrix<T,0,0>& a, const SimpleMatrix<T,0,0>& b) {
    BOOST_CHECK_EQUAL(a.rows(), b.rows());
    BOOST_CHECK_EQUAL(a.cols(), b.cols());
    
    T residual = 0;
    for (index_t r = 0; r < a.rows(); ++r) {
        for (index_t c = 0; c < a.cols(); ++c) {
            T z = a[r][c] - b[r][c];
            residual += z * z;
        }
    }
    return std::sqrt(residual);
}

template <typename T>
void run_cholesky(index_t sz, rng_t* rng) {
    SimpleMatrix<T,0,0> mx = random_matrix<T>(sz, rng);
    SimpleMatrix<T,0,0> mxT(sz, sz);
    SimpleMatrix<T,0,0> A(sz, sz);
    SimpleMatrix<T,0,0> C(sz, sz);
    
    // make mx a positive definite matrix `A` by taking mx * mx^T:
    transpose(&mxT, mx);
    mul(&A, mx, mxT);
    
    // cholesky decompose `A` into `mx`.
    mtxcopy(&mx, A);
    BOOST_CHECK(cholesky(&mx));
    
    // confirm `mx^T * mx = A`
    transpose(&mxT, mx);
    mul(&C, mx, mxT);
    T rms = matrix_diff(C, A);
    BOOST_CHECK_SMALL(rms, (T)1e-5);
}


BOOST_AUTO_TEST_SUITE(cholesky)


BOOST_AUTO_TEST_CASE(verify_cholesky) {
    rng_t rng(11937294775LL);
    std::uniform_int_distribution<> rnd_int(2, 16);
    const index_t n = 50000;
    for (index_t i = 0; i < n; ++i) {
        index_t k = rnd_int(rng);
        run_cholesky<double>(k, &rng);
    }
}


BOOST_AUTO_TEST_SUITE_END()
