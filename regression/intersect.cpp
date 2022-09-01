#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Intersect

#include <random>
#include <boost/test/unit_test.hpp>

#include <geomc/shape/Intersect.h>

#include "shape_generation.h"

template <typename T, index_t N>
void random_box(OrientedRect<T,N>* r, rng_t* rng) {
    std::uniform_real_distribution<T> range{2,8};
    
    rnd(rng, &r->xf);
    r->xf *= scale(Vec<T,N>(range(*rng)));
    r->xf *= translation(rnd<T,N>(rng).unit());
    ShapeSampler<Rect<T,N>> box{{-1,1}};
    r->shape = Rect<T,N>::spanning_corners(box(rng), box(rng));
}

template <typename T, index_t N>
void get_box_corners(const OrientedRect<T,N>& r, Vec<T,N> p[1 << N]) {
    Vec<T,N> extreme[2] = { r.shape.lo, r.shape.hi };
    
    for (index_t i = 0; i < (1 << N); ++i) {
        Vec<T,N> v;
        for (index_t axis = 0; axis < N; ++axis) {
            v[axis] = extreme[((i >> axis) & 1)][axis];
        }
        p[i] = r.xf * v;
    }
}

// this is tautological for N > 3. OBBs use GJK directly!
template <typename T, index_t N> void test_box_gjk_matches_sat(rng_t* rng, index_t iters) {
    const index_t n = (index_t)std::ceil(std::sqrt(iters));
    OrientedRect<T,N>* boxes = new OrientedRect<T,N>[n];
    for (index_t i = 0; i < n; i++) {
        random_box(boxes + i, rng);
    }
    constexpr index_t n_corners = 1 << N;
    
    Vec<T,N> d;
    index_t i0 = 0;
    index_t i1 = 0;
    index_t failures = 0;
    index_t positive = 0;
    index_t negative = 0;
    for (index_t j = 0; j < iters; j++) {
        i0 = (i0 + 1) % n;
        if (i0 == 0) i1 = (i1 + 1) % n;
        Vec<T,N> b0[n_corners];
        Vec<T,N> b1[n_corners];
        get_box_corners(boxes[i0], b0);
        get_box_corners(boxes[i1], b1);
        bool gjk = gjk_intersect(b0, n_corners, b1, n_corners, &d);
        bool SAT = boxes[i0].intersects(boxes[i1]);
        if (gjk != SAT) {
            failures++;
        }
        if (SAT) positive++;
        else negative++;
    }
    
    delete [] boxes;
    std::cout << "tested " << iters << " shapes. " << positive << " intersecting, ";
    std::cout << negative << " disjoint" << std::endl;
    BOOST_CHECK_EQUAL(failures, 0);
}

BOOST_AUTO_TEST_SUITE(intersect)

BOOST_AUTO_TEST_CASE(test_gjk_sat_agree) {
    const index_t N_TESTS = 1000000;
    test_box_gjk_matches_sat<float, 2>(&rng, N_TESTS);
    test_box_gjk_matches_sat<float, 3>(&rng, N_TESTS);
    test_box_gjk_matches_sat<double, 2>(&rng, N_TESTS);
    test_box_gjk_matches_sat<double, 3>(&rng, N_TESTS);
}

BOOST_AUTO_TEST_SUITE_END()