#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE LinearSolve

#include <iostream>
#include <random>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <geomc/linalg/LUDecomp.h>
#include <geomc/linalg/Orthogonal.h>

using namespace geom;
using namespace std;

typedef std::mt19937_64 rng_t;

// test: check that PLU = M, actually
//       (is P what it says it is?)


template <typename T>
void randomize(T* v, index_t n, rng_t* rng) {
    typedef std::normal_distribution<T> d_normal_t;
    d_normal_t gauss = d_normal_t(0., 1); // (ctr, variance)
    
    for (index_t i = 0; i < n; ++i) {
        v[i] = gauss(*rng);
    }
}


template <typename T, index_t N>
void solve_vectors(rng_t* rng) {
    Vec<T,N> bases[N]; // destroyed by the factorization
    Vec<T,N> vi[N];
    Vec<T,N> b;
    Vec<T,N> x;
    randomize<T>(bases[0].begin(), N * N, rng);
    randomize<T>(b.begin(), N, rng);
    std::copy(bases, bases + N, vi); // copy the bases, so we can check them
    
    // we should be able to write b as x[i] * bases[i].
    if (linear_solve(bases, &x, b, 0)) {
        Vec<T,N> b0;
        for (index_t i = 0; i < N; ++i) {
            b0 += x[i] * vi[i];
        }
        // check that the components of b0 and b are close.
        for (index_t i = 0; i < N; ++i) {
            BOOST_CHECK_CLOSE(b[i], b0[i], 1e-5);
        }
    } else {
        // ↓ very unlikely
        std::cerr << "Singular test matrix\n";
    }
}


template <typename T, index_t N>
void exercise_solve_vectors(rng_t* rng, index_t trials) {
    for (index_t i = 0; i < trials; ++i) {
        solve_vectors<T,N>(rng);
    }
}


template <typename T>
void exercise_plu(rng_t* rng, index_t n, index_t trials) {
    SimpleMatrix<T,0,0> mx(n, n);
    SimpleMatrix<T,0,0>  L(n,n);
    SimpleMatrix<T,0,0>  U(n,n);
    SimpleMatrix<T,0,0> LU(n,n);
    PLUDecomposition<T,0,0> plu(mx);
    
    for (index_t i = 0; i < trials; ++i) {
        randomize(mx.begin(), n * n, rng);
        plu.decompose(mx);
        plu.get_L(&L);
        plu.get_U(&U);
        mul(&LU, L, U);     // LU <- L * U
        mul(&L, plu.P, mx); // L  <- P * M
        // check that LU = PM
        for (index_t j = 0; j < n * n; ++j) {
            BOOST_CHECK_CLOSE(LU.begin()[j], L.begin()[j], 1e-5);
        }
    }
}


template <typename T, index_t N>
void exercise_orthogonalize(rng_t* rng, index_t trials) {
    Vec<T,N> vs[N];
    for (index_t i = 0; i < trials; ++i) {
        randomize(vs[0].begin(), N * N, rng);
        orthogonalize(vs, N);
        // check that each of the orthogonalized bases have ~= 0 dot
        // product with all the others:
        for (index_t j = 0; j < N - 1; ++j) {
            for (index_t k = j + 1; k < N; ++k) {
                BOOST_CHECK_SMALL(vs[j].dot(vs[k]), 1e-5);
            }
        }
    }
}


template <typename T, index_t N>
void exercise_nullspace(rng_t* rng, index_t trials) {
    auto urnd = std::uniform_int_distribution<index_t>(1, N - 1);
    Vec<T,N> bases[N];
    for (index_t i = 0; i < trials; ++i) {
        index_t n = urnd(*rng);
        randomize(bases[0].begin(), n * N, rng);
        nullspace(bases, n, bases + n);
        
        // verify that each of the null bases have
        // zero dot product with the source bases.
        // (they are not guaranteed to be orthogonal to each other)
        for (index_t j = 0; j < n; ++j) {
            for (index_t k = n; k < N; ++k) {
                BOOST_CHECK_SMALL(bases[j].dot(bases[k]), 1e-5);
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE(linear_solve_tests)


BOOST_AUTO_TEST_CASE(verify_nullspace) {
    rng_t rng(16512485420001724907ULL);
    exercise_nullspace<double,  2>(&rng, 250);
    exercise_nullspace<double,  3>(&rng, 250);
    exercise_nullspace<double,  4>(&rng, 250);
    exercise_nullspace<double,  5>(&rng, 250);
    exercise_nullspace<double, 10>(&rng, 250);
}


BOOST_AUTO_TEST_CASE(verify_orthogonal) {
    rng_t rng(15794404771588593305ULL);
    exercise_orthogonalize<double,  2>(&rng, 250);
    exercise_orthogonalize<double,  3>(&rng, 250);
    exercise_orthogonalize<double,  4>(&rng, 250);
    exercise_orthogonalize<double,  5>(&rng, 250);
    exercise_orthogonalize<double, 10>(&rng, 250);
}


BOOST_AUTO_TEST_CASE(verify_PLU) {
    rng_t rng(1013126187766094264ULL);
    exercise_plu<double>(&rng, 2,  250);
    exercise_plu<double>(&rng, 3,  250);
    exercise_plu<double>(&rng, 4,  250);
    exercise_plu<double>(&rng, 5,  250);
    exercise_plu<double>(&rng, 7,  250);
    exercise_plu<double>(&rng, 10, 100);
}


BOOST_AUTO_TEST_CASE(linear_solve_tests) {
    rng_t rng(7301667549950575693ULL);
    exercise_solve_vectors<double,2>(&rng, 250);
    exercise_solve_vectors<double,3>(&rng, 250);
    exercise_solve_vectors<double,4>(&rng, 250);
    exercise_solve_vectors<double,5>(&rng, 250);
    exercise_solve_vectors<double,6>(&rng, 250);
    exercise_solve_vectors<double,7>(&rng, 250);
}


BOOST_AUTO_TEST_SUITE_END()

