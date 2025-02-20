#define TEST_MODULE_NAME Perlin

#include <pcg_random.hpp>
#include <gtest/gtest.h>
#include <geomc/shape/Sphere.h>
#include <geomc/random/SampleGeometry.h>
#include <geomc/function/PerlinNoise.h>

using namespace geom;
using namespace std;

using rng_t = pcg64;


// todo: should figure out a way to get a copy of a PerlinNoise in Dual form, 
//       and verify the gradient that way.


// return rms(error), max(error)
template <typename T, index_t N>
std::pair<T,T> perlin_gradient(rng_t& rng, index_t n_trials) {
    PerlinNoise<T,N> pn(rng);
    const T eps = 0.00001;
    
    T err_sq  = (T)0;
    T err_max = (T)0;
    
    SampleShape<Sphere<T,N>> smp_sphere {};
    
    for (index_t i = 0; i < n_trials; ++i) {
        Vec<T,N>  x = 512 * smp_sphere(rng); // :G
        auto   x_dx = pn.gradient(x);
        Vec<T,N> g;
        
        // finite difference the gradient
        for (index_t axis = 0; axis < N; axis++) {
            Vec<T,N> dx;
            dx[axis] = eps;
            g[axis]  = (pn.eval(x + dx) - pn.eval(x - dx)) / (2 * eps);
        }
        
        // gradient()'s opinion on f(x) should be the same as eval()'s.
        EXPECT_NEAR(x_dx.first, pn.eval(x), eps);
        
        g      -= x_dx.second;
        T e     = g.dot(g);
        err_sq += e;
        err_max = std::max(e, err_max);
    }
    
    return std::pair<T,T>(
        std::sqrt(err_sq) / n_trials,
        std::sqrt(err_max));
}


TEST(TEST_MODULE_NAME, test_perlin_gradient) {
    rng_t rng {1017381749271967481LL};
    std::pair<double, double> k;
    k = perlin_gradient<double,2>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 5e-11);
    EXPECT_NEAR(k.second, 0, 1e-8);
    k = perlin_gradient<double,3>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 5e-11);
    EXPECT_NEAR(k.second, 0, 1e-8);
    k = perlin_gradient<double,4>(rng, 10000);
    EXPECT_NEAR(k.first,  0, 5e-11);
    EXPECT_NEAR(k.second, 0, 1e-8);
}
