#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Perlin

#include <boost/test/unit_test.hpp>
#include <geomc/function/PerlinNoise.h>
#include <geomc/random/MTRand.h>
#include <geomc/random/RandomTools.h>

using namespace geom;
using namespace std;


// todo: should figure out a way to get a copy of a PerlinNoise in Dual form, 
//       and verify the gradient that way.


// return rms(error), max(error)
template <typename T, index_t N>
std::pair<T,T> perlin_gradient(Random* rng, index_t n_trials) {
    PerlinNoise<T,N> pn(rng);
    Sampler<T> smp(rng);
    const T eps = 0.00001;
    
    T err_sq  = (T)0;
    T err_max = (T)0;
    
    for (index_t i = 0; i < n_trials; ++i) {
        Vec<T,N>  x = 512 * smp.template solidball<N>(); // :G
        auto   x_dx = pn.gradient(x);
        Vec<T,N> g;
        
        // finite difference the gradient
        for (index_t axis = 0; axis < N; axis++) {
            Vec<T,N> dx;
            dx[axis] = eps;
            g[axis]  = (pn.eval(x + dx) - pn.eval(x - dx)) / (2 * eps);
        }
        
        // gradient()'s opinion on f(x) should be the same as eval()'s.
        BOOST_CHECK_CLOSE(x_dx.first, pn.eval(x), eps);
        
        g      -= x_dx.second;
        T e     = g.dot(g);
        err_sq += e;
        err_max = std::max(e, err_max);
    }
    
    return std::pair<T,T>(
        std::sqrt(err_sq) / n_trials,
        std::sqrt(err_max));
}


BOOST_AUTO_TEST_SUITE(perlin_noise)


BOOST_AUTO_TEST_CASE(test_perlin_gradient) {
    MTRand rng = MTRand(1017381749271967481LL);
    std::pair<double, double> k;
    k = perlin_gradient<double,2>(&rng, 10000);
    BOOST_CHECK_SMALL(k.first,  5e-11);
    BOOST_CHECK_SMALL(k.second, 1e-8);
    k = perlin_gradient<double,3>(&rng, 10000);
    BOOST_CHECK_SMALL(k.first,  5e-11);
    BOOST_CHECK_SMALL(k.second, 1e-8);
    k = perlin_gradient<double,4>(&rng, 10000);
    BOOST_CHECK_SMALL(k.first,  5e-11);
    BOOST_CHECK_SMALL(k.second, 1e-8);
}


BOOST_AUTO_TEST_SUITE_END()
