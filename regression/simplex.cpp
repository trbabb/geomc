#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Simplex

#include <iostream>
#include <random>
#include <pcg_random.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <geomc/shape/Simplex.h>
#include <geomc/linalg/Matrix.h>

#include "shape_generation.h"

// increase for good coverage.
// set to 1 for debugging.
#define N_TESTS 100'000

using namespace geom;
using namespace std;

typedef pcg64 rng_t;

template <typename T>
void randomize(T* v, index_t n, rng_t* rng) {
    typedef std::normal_distribution<T> d_normal_t;
    d_normal_t gauss = d_normal_t(0., 1); // (ctr, variance)
    
    for (index_t i = 0; i < n; ++i) {
        v[i] = gauss(*rng);
    }
}

template <typename T, index_t N>
bool validate_simplex_measure(const Simplex<T,N>& s) {
    index_t n = s.n - 1;
    T vol = s.measure();
    int k = 1;
    Vec<T,N> bases[N];
    for (index_t i = 0; i < n; ++i) {
        bases[i] = s[i + 1] - s[0];
        k *= (i + 1);
    }
    // volume is sqrt of the determinant of the gram matrix
    // divided by n factorial.
    T b[N * N];
    WrapperMatrix<T,0,0> m  {b, n, n};
    WrapperMatrix<T,0,N,MatrixLayout::ROW_MAJOR> m0 {bases[0].begin(), n};
    WrapperMatrix<T,N,0,MatrixLayout::COL_MAJOR> m1 {bases[0].begin(), n};
    // (n x n) = (n x N) * (N x n)
    mul(&m, m0, m1);
    T vol2 = std::sqrt(std::abs(det_destructive(b, n))) / k;
    
    // std::cout << "error: " << std::abs(vol2 / vol - 1) << "\n";
    
    // BOOST_CHECK_SMALL(vol2 / vol - 1, 1e-4);
    // if (std::abs(vol2 / vol - 1) > 1e-5) {
    //     // std::cout << "incorrect | N: " << N << " n: " << n << "\n";
    //     // std::cerr << "  vol: " << vol << " vol2: " << vol2 << "\n";
    //     // std::cerr << "    N: " << N   << "    n: " << n    << "\n";
    // } else {
    //     std::cout << "correct | N: " << N << " n: " << n << "\n";
    // }
    return std::abs(vol2 / vol - 1) < 1e-5;
}


template <typename T, index_t N>
void exercise_simplex_measure(rng_t* rng, index_t trials) {
    auto urnd = std::uniform_int_distribution<index_t>(1, N + 1);
    Simplex<T,N> s;
    int ok[N + 1] = {0};
    for (index_t i = 0; i < trials; ++i) {
        index_t n = urnd(*rng);
        randomize(s.pts[0].begin(), n * N, rng);
        s.n = n;
        ok[n] += validate_simplex_measure(s);
    }
    std::cout << "--------------------- " << N << "D ---------------------\n";
    for (index_t i = 1; i <= N; ++i) {
        std::cout << "n = " << i << " | " << ok[i] / (float) trials << " success\n";
    }
    // std::cout << "---------> N = " << N << " | " << ok / (float) trials << " success\n";
}


template <typename T, index_t N>
void test_simplex_projection(rng_t* rng, const Simplex<T,N>& s, index_t trials) {
    auto bb   = s.bounds();
    auto dims = bb.dimensions();
    auto ctr  = bb.center();
    // no degenerate boxes pls
    for (index_t i = 0; i < N; ++i) [[unlikely]] {
        if (dims[i] == 0) dims[i] = 1;
    }
    for (index_t i = 0; i < trials; ++i) {
        Simplex<T,N> ss_proj;
        Simplex<T,N> ss_clip;
        Vec<T,N> p  = rnd<T,N>(rng) * dims + ctr;
        Vec<T,N> pp = s.project(p, &ss_proj);
        Vec<T,N> pc = s.clip(p, &ss_clip);
        bool inside = s.contains(p);
        if (inside) {
            // `p` is inside a full-volume simplex.
            // no points should have been excluded;
            // the "projected-to" simplex should be the same as the original
            BOOST_CHECK(ss_clip == s);
            // the target simplex is the same as the original
            BOOST_CHECK_EQUAL(ss_clip.n, s.n);
            // the clipped point should be precisely the original point
            BOOST_CHECK_EQUAL(p, pc);
            // the projected point should be an n-1 dimensional face
            BOOST_CHECK_EQUAL(ss_proj.n, N);
        } else {
            // if the point is not inside the simplex, then projected point
            // and the clipped point are the same
            BOOST_CHECK_EQUAL(pp, pc);
            // the direction from the surface pt to the original point
            // should be orthogonal to the simplex face
            Vec<T,N> v = p - pp;
            for (index_t j = 1; j < ss_clip.n; ++j) {
                Vec<T,N> b = ss_clip.pts[j] - ss_clip.pts[0];
                BOOST_CHECK_SMALL(b.dot(v), 1e-5);
            }
            bool facing = ss_clip.is_facing(p);
            BOOST_CHECK(facing);
        }
    }
}


template <typename T, index_t N>
void exercise_simplex_projection(rng_t* rng, index_t trials) {
    for (index_t i = 0; i < trials; ++i) {
        Simplex<T,N> s = RandomShape<Simplex<T,N>>::rnd_shape(rng);
        std::uniform_int_distribution<> d(1,N+1);
        s.n = d(*rng); // make the simplex have a random number of verts
        test_simplex_projection(rng, s, 10);
    }
}


BOOST_AUTO_TEST_SUITE(simplex_tests)

/*
BOOST_AUTO_TEST_CASE(verify_simplex_measure) {
    rng_t rng(16512485420001724907ULL);
    exercise_simplex_measure<double,  2>(&rng, N_TESTS);
    exercise_simplex_measure<double,  3>(&rng, N_TESTS);
    exercise_simplex_measure<double,  4>(&rng, N_TESTS);
    exercise_simplex_measure<double,  5>(&rng, N_TESTS);
    exercise_simplex_measure<double, 10>(&rng, N_TESTS);
}
*/


BOOST_AUTO_TEST_CASE(simplex_projection) {
    rng_t rng(2532266933125789701ULL);
    std::cout << "computing 2D..." << std::endl;
    exercise_simplex_projection<double,2>(&rng, N_TESTS);
    std::cout << "computing 3D..." << std::endl;
    exercise_simplex_projection<double,3>(&rng, N_TESTS);
    std::cout << "computing 4D..." << std::endl;
    exercise_simplex_projection<double,4>(&rng, N_TESTS);
    std::cout << "computing 5D..." << std::endl;
    exercise_simplex_projection<double,5>(&rng, std::max(N_TESTS/10,  1));
    std::cout << "computing 7D..." << std::endl;
    exercise_simplex_projection<double,7>(&rng, std::max(N_TESTS/100, 1));
}


BOOST_AUTO_TEST_SUITE_END()

