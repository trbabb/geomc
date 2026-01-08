#define TEST_MODULE_NAME Similarity

#include <random>
#include <pcg_random.hpp>
#include <gtest/gtest.h>
#include <geomc/linalg/Similarity.h>
#include <geomc/linalg/AffineTransform.h>
#include <geomc/random/SampleVector.h>

#define N_TESTS 65'535

using namespace geom;
using namespace std;

typedef pcg64 rng_t;

template <typename T, index_t N>
Similarity<T,N> random_similarity(rng_t* rng);

// 2D specialization
template <>
Similarity<double,2> random_similarity<double,2>(rng_t* rng) {
    std::normal_distribution<double> gauss(0, 1);
    std::uniform_real_distribution<double> angle_dist(0, 2 * std::numbers::pi);

    double sx = gauss(*rng);
    Rotation<double,2> rx(angle_dist(*rng));
    Vec<double,2> tx = random_gaussian<double,2>(*rng);

    return Similarity<double,2>(sx, rx, tx);
}

// 3D specialization
template <>
Similarity<double,3> random_similarity<double,3>(rng_t* rng) {
    std::normal_distribution<double> gauss(0, 1);
    std::uniform_real_distribution<double> angle_dist(0, 2 * std::numbers::pi);

    double sx = gauss(*rng);
    Quat<double> q = random_unit<double,4>(*rng);
    Rotation<double,3> rx(q);
    Vec<double,3> tx = random_gaussian<double,3>(*rng);

    return Similarity<double,3>(sx, rx, tx);
}

// Test that the AffineTransform conversion produces equivalent transformations
template <typename T, index_t N>
void test_affine_conversion(rng_t* rng) {
    constexpr index_t N_VECS = 100;

    Similarity<T,N> sim = random_similarity<T,N>(rng);
    AffineTransform<T,N> xf = sim;

    // Test zero vector
    Vec<T,N> zero;
    Vec<T,N> sim_zero = sim * zero;
    Vec<T,N> xf_zero = xf * zero;
    for (index_t i = 0; i < N; ++i) {
        EXPECT_NEAR(sim_zero[i], xf_zero[i], 1e-5);
    }

    // Test random vectors
    for (index_t j = 0; j < N_VECS; ++j) {
        Vec<T,N> v = random_gaussian<T,N>(*rng);
        Vec<T,N> sim_v = sim * v;
        Vec<T,N> xf_v = xf * v;

        for (index_t i = 0; i < N; ++i) {
            EXPECT_NEAR(sim_v[i], xf_v[i], 1e-5);
        }
    }
}

template <typename T, index_t N>
void exercise_affine_conversion(rng_t* rng, index_t trials) {
    for (index_t i = 0; i < trials; ++i) {
        test_affine_conversion<T,N>(rng);
    }
}

// Test that composition of Similarities is equivalent to composition of their AffineTransforms
template <typename T, index_t N>
void test_composition(rng_t* rng, index_t n_compose) {
    constexpr index_t N_VECS = 50;

    // Generate n_compose random similarities
    Similarity<T,N> sims[5];

    for (index_t i = 0; i < n_compose; ++i) {
        sims[i] = random_similarity<T,N>(rng);
    }

    // Compose them
    Similarity<T,N> sim_composed;
    AffineTransform<T,N> xf_composed;

    for (index_t i = 0; i < n_compose; ++i) {
        AffineTransform<T,N> xf_next = sims[i];
        sim_composed = sim_composed * sims[i];
        xf_composed = xf_composed * xf_next;
    }

    // Test that they produce the same results on random vectors
    for (index_t j = 0; j < N_VECS; ++j) {
        Vec<T,N> v = random_gaussian<T,N>(*rng);
        Vec<T,N> sim_v = sim_composed * v;
        Vec<T,N> xf_v = xf_composed * v;

        for (index_t i = 0; i < N; ++i) {
            EXPECT_NEAR(sim_v[i], xf_v[i], 1e-5);
        }
    }

    // Also test the zero vector
    Vec<T,N> zero;
    Vec<T,N> sim_zero = sim_composed * zero;
    Vec<T,N> xf_zero = xf_composed * zero;
    for (index_t i = 0; i < N; ++i) {
        EXPECT_NEAR(sim_zero[i], xf_zero[i], 1e-5);
    }
}

template <typename T, index_t N>
void exercise_composition(rng_t* rng, index_t trials) {
    std::uniform_int_distribution<index_t> n_compose_dist(2, 5);

    for (index_t i = 0; i < trials; ++i) {
        index_t n_compose = n_compose_dist(*rng);
        test_composition<T,N>(rng, n_compose);
    }
}

TEST(TEST_MODULE_NAME, verify_affine_conversion_2d) {
    rng_t rng(0xee5ecfb650cdccdc);
    exercise_affine_conversion<double,2>(&rng, N_TESTS);
}

TEST(TEST_MODULE_NAME, verify_affine_conversion_3d) {
    rng_t rng(0x17c73a76e09ad68);
    exercise_affine_conversion<double,3>(&rng, N_TESTS);
}

TEST(TEST_MODULE_NAME, verify_composition_2d) {
    rng_t rng(0x9c503a5574aecf45);
    exercise_composition<double,2>(&rng, N_TESTS);
}

TEST(TEST_MODULE_NAME, verify_composition_3d) {
    rng_t rng(0x9c503a5574aecf45);
    exercise_composition<double,3>(&rng, N_TESTS);
}
