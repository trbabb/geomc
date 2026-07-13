#define TEST_MODULE_NAME SampleGeometry

#include <cmath>

#include <pcg_random.hpp>
#include <gtest/gtest.h>

#include <geomc/random/SampleGeometry.h>
#include <geomc/shape/Sphere.h>
#include <geomc/shape/Hollow.h>

using namespace geom;

using rng_t = pcg64;

// Points sampled from the surface of a sphere (Hollow<Sphere>) must lie exactly
// on the sphere: |sample - center| == radius. This exercises the surface sampler
// across both the low-dimensional rejection path (N <= 3) and the gaussian path
// (N >= 4).
template <typename T, index_t N>
void check_sphere_surface(rng_t& rng, index_t n_trials) {
    Sphere<T,N> sphere {VecType<T,N>((T)3), (T)2}; // center (3,3,...), radius 2
    SampleShape<Hollow<Sphere<T,N>>> sampler {Hollow<Sphere<T,N>>(sphere)};
    const T eps = 100 * std::numeric_limits<T>::epsilon();
    for (index_t i = 0; i < n_trials; ++i) {
        VecType<T,N> p = sampler(rng);
        T r = std::sqrt(PointType<T,N>::mag2(p - sphere.center));
        EXPECT_NEAR(r, sphere.radius, eps * sphere.radius);
    }
}

TEST(TEST_MODULE_NAME, sphere_surface_is_unit) {
    rng_t rng {8675309ULL};
    check_sphere_surface<double,2>(rng, 4096);
    check_sphere_surface<double,3>(rng, 4096);
    check_sphere_surface<double,4>(rng, 4096); // gaussian path
    check_sphere_surface<double,5>(rng, 4096);
    check_sphere_surface<float,2>(rng, 4096);
    check_sphere_surface<float,3>(rng, 4096);
    check_sphere_surface<float,4>(rng, 4096);
}
