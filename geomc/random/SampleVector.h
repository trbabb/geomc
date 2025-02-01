#pragma once

#include <random>
#include <geomc/linalg/Vec.h>

namespace geom {

/**
 * @addtogroup random
 * @{
 */

/// A default random number generator type.
using DefaultLCG = std::linear_congruential_engine<
    // Knuth's MMIX LCG parameters
    uint64_t,
    6364136223846793005ULL,
    1442695040888963407ULL,
    0
>;

/// Create a new random number generator with a nondeterministic seed.
inline DefaultLCG create_rng() {
    using T = DefaultLCG::result_type;
    using Seeder = std::uniform_int_distribution<T>;
    std::random_device rd;
    return DefaultLCG(Seeder{}(rd));
}

/**
 * @brief Generate a random vector drawn from a multivariate gaussian distribution
 * with mean 0 and variance 1.
 */
template <typename T, index_t N, typename Generator>
inline Vec<T,N> random_gaussian(Generator& rng) {
    std::normal_distribution<T> gauss(0, 1);
    Vec<T,N> v;
    for (index_t i = 0; i < N; ++i) {
        v[i] = gauss(rng);
    }
    return v;
}

/**
 * @brief Generate a random vector with unit length.
 */
template <typename T, index_t N, typename Generator>
inline Vec<T,N> random_unit(Generator& rng) {
    Vec<T,N> p;
    if constexpr (N <= 3) {
        // rejection sampling
        std::uniform_real_distribution<T> unif {-1,1};
        do {
            // generate a point in the signed unit box
            for (index_t i = 0; i < N; ++i) {
                p[i] = unif(rng);
            }
            // if the point is outside the unit sphere, reject it and try again
        } while (p.mag2() > 1);
        return p.unit();
    } else {
        // draw a multivariate gaussian, then project it onto the sphere
        std::normal_distribution<T> gauss(0, 1);
        for (index_t i = 0; i < N; ++i) {
            p[i] = gauss(rng);
        }
    }
    return p.unit();
}

/// @}

}  // namespace geom
