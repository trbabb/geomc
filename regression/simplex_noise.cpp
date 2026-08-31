#define TEST_MODULE_NAME SimplexNoise

#include <pcg_random.hpp>
#include <gtest/gtest.h>
#include <geomc/shape/Sphere.h>
#include <geomc/random/SampleGeometry.h>
#include <geomc/function/SimplexNoise.h>

using namespace geom;
using namespace std;

using rng_t = pcg64;


// return rms(error), max(error) between the analytic gradient and a central
// finite difference of the noise value. mirrors regression/perlin.cpp.
template <typename T, index_t N>
std::pair<T,T> simplex_gradient(rng_t& rng, index_t n_trials) {
    SimplexNoise<T,N> sn(rng);
    const T eps = 0.00001;

    T err_sq  = (T)0;
    T err_max = (T)0;

    SampleShape<Sphere<T,N>> smp_sphere {};

    for (index_t i = 0; i < n_trials; ++i) {
        Vec<T,N>  x = 512 * smp_sphere(rng); // :G
        auto   x_dx = sn.gradient(x);
        Vec<T,N> g;

        // finite difference the gradient
        for (index_t axis = 0; axis < N; axis++) {
            Vec<T,N> dx;
            dx[axis] = eps;
            g[axis]  = (sn.eval(x + dx) - sn.eval(x - dx)) / (2 * eps);
        }

        // gradient()'s opinion on f(x) should be the same as eval()'s.
        EXPECT_NEAR(x_dx.first, sn.eval(x), eps);

        g      -= x_dx.second;
        T e     = g.dot(g);
        err_sq += e;
        err_max = std::max(e, err_max);
    }

    return std::pair<T,T>(
        std::sqrt(err_sq) / n_trials,
        std::sqrt(err_max));
}


// Walk a dense line and confirm the gradient never jumps: a C1 field has a
// bounded finite-difference of the analytic gradient. A merely-C0 (faceted)
// field would show unbounded spikes at simplex boundaries.
template <typename T, index_t N>
T max_gradient_step(rng_t& rng, index_t steps) {
    SimplexNoise<T,N> sn(rng);
    Vec<T,N> o, dir;
    for (index_t k = 0; k < N; ++k) { o[k] = 20 * (T)(k + 1); dir[k] = (T)(k + 1); }
    dir = dir / std::sqrt(dir.mag2());
    const T len = 60;
    const T h   = len / steps;
    Vec<T,N> prev;
    T max_dg = 0;
    for (index_t i = 0; i < steps; ++i) {
        Vec<T,N> x = o + dir * (h * i);
        Vec<T,N> g = sn.gradient(x).second;
        if (i > 0) max_dg = std::max(max_dg, std::sqrt((g - prev).mag2()) / h);
        prev = g;
    }
    return max_dg;
}


TEST(TEST_MODULE_NAME, gradient_matches_finite_difference) {
    rng_t rng {1017381749271967481LL};
    std::pair<double, double> k;
    k = simplex_gradient<double,2>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 1e-9);
    EXPECT_NEAR(k.second, 0, 1e-6);
    k = simplex_gradient<double,3>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 1e-9);
    EXPECT_NEAR(k.second, 0, 1e-6);
    k = simplex_gradient<double,4>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 1e-9);
    EXPECT_NEAR(k.second, 0, 1e-6);
    k = simplex_gradient<double,5>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 1e-9);
    EXPECT_NEAR(k.second, 0, 1e-6);
}


TEST(TEST_MODULE_NAME, gradient_is_continuous) {
    rng_t rng {8899223344556677LL};
    // C1 continuity: the gradient's rate of change along a line stays bounded.
    // (A faceted, C0 field would produce spikes orders of magnitude larger.)
    EXPECT_LT((max_gradient_step<double,2>(rng, 200000)), 100.0);
    EXPECT_LT((max_gradient_step<double,3>(rng, 200000)), 100.0);
    EXPECT_LT((max_gradient_step<double,4>(rng, 200000)), 100.0);
}
