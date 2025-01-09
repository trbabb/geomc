#define TEST_MODULE_NAME Intersect

// set this to 1 to get diagnostic text file dumps of failure cases
#define DEBUG_INTERSECTION 0

#include <random>
#include <gtest/gtest.h>

#include <geomc/shape/Intersect.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/Transformed.h>

#include "shape_generation.h"

#if DEBUG_INTERSECTION
constexpr bool DO_DEBUG = true;
#else
constexpr bool DO_DEBUG = false;
#endif
    
template <typename T, index_t N>
struct PointHull : public Convex<T,N,PointHull<T,N>> {
    const Vec<T,N>* pts;
    index_t n;
    
    PointHull(const Vec<T,N> *pts, index_t n):pts(pts),n(n) {}
    
    Vec<T,N> convex_support(Vec<T,N> d) const {
        T largest_dot = std::numeric_limits<T>::lowest();
        Vec<T,N> pt;
        for (index_t i = 0; i < n; i++) {
            T a = d.dot(pts[i]);
            if (a > largest_dot) {
                largest_dot = a;
                pt = pts[i];
            }
        }
        return pt;
    }
};

template <typename T, index_t N>
void random_box(AffineBox<T,N>* r, rng_t* rng) {
    std::uniform_real_distribution<T> range{2,8};
    
    rnd(rng, &r->xf);
    r->xf *= scale(Vec<T,N>(range(*rng)));
    r->xf *= translation(rnd<T,N>(rng).unit());
    ShapeSampler<Rect<T,N>> box{{-1,1}};
    r->shape = Rect<T,N>::spanning_corners(box(rng), box(rng));
}

template <typename T, index_t N>
void get_box_corners(const AffineBox<T,N>& r, Vec<T,N> p[1 << N]) {
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
template <typename T, index_t N>
void test_box_gjk_matches_sat(rng_t* rng, index_t iters) {
    const index_t n = (index_t)std::ceil(std::sqrt(iters));
    AffineBox<T,N>* boxes = new AffineBox<T,N>[n];
    for (index_t i = 0; i < n; i++) {
        random_box(boxes + i, rng);
    }
    constexpr index_t n_corners = 1 << N;
    
    std::vector<index_t> iter_hist;
    
    char fname[] = "gjk_D_.txt";
    fname[3] = '0' + N;
    fname[5] = std::is_same<T,float>::value ? 'f' : 'd';
    FILE* f = DO_DEBUG ? fopen(fname, "w") : nullptr;
    
    index_t i0 = 0;
    index_t i1 = 0;
    index_t failures = 0;
    index_t positive = 0;
    index_t negative = 0;
    index_t pos_fails = 0;
    index_t neg_fails = 0;
    index_t pos_results = 0;
    index_t neg_results = 0;
    index_t ident_failures = 0;
    index_t degen_failures = 0;
    index_t max_iters = 0;
    for (index_t j = 0; j < iters; j++) {
        i0 = (i0 + 1) % n;
        if (i0 == 0) i1 = (i1 + 1) % n;
        Vec<T,N> b0[n_corners];
        Vec<T,N> b1[n_corners];
        get_box_corners(boxes[i0], b0);
        get_box_corners(boxes[i1], b1);
        Intersector<T,N> intersector;
        bool gjk = intersector.intersects(
            as_any_convex(PointHull<T,N>(b0, n_corners)), 
            as_any_convex(PointHull<T,N>(b1, n_corners))
        );
        bool SAT = boxes[i0].intersects(boxes[i1]);
        index_t iters = intersector.iterations;
        if (iter_hist.size() < iters + 1) iter_hist.resize(iters + 1);
        iter_hist[iters] += 1;
        max_iters = std::max(max_iters, iters);
        if (gjk != SAT) {
            failures++;
            if (SAT) pos_fails++;
            else neg_fails++;
            if (i0 == i1) ident_failures++;
            if (intersector.was_degenerate) degen_failures++;
            if constexpr (DO_DEBUG) {
                if (failures < 30) {
                    // emit a limited number of failure cases
                    emit_hull(f, b0, n_corners);
                    emit_hull(f, b1, n_corners);
                    intersector = {};
                    intersector.debug_file = f;
                    intersector.intersects(
                        as_any_convex(PointHull<T,N>(b0, n_corners)), 
                        as_any_convex(PointHull<T,N>(b1, n_corners))
                    );
                    fprintf(f, "========\n");
                }
            }
        }
        if (gjk) pos_results++;
        else neg_results++;
        
        if (SAT) positive++;
        else negative++;
    }
    
    if (f) fclose(f);
    
    delete [] boxes;
    std::cout << "tested " << iters << " " << N << "D" << typeid(T).name() << " shapes. ";
    std::cout << positive << " intersecting, " << negative << " disjoint" << std::endl;
    if (failures != 0) {
        std::cout << "  actually-intersecting failures: "  << pos_fails << std::endl;
        std::cout << "  actually-disjoint failures: "      << neg_fails << std::endl;
        std::cout << "  identical-object failures: "       << ident_failures << std::endl;
        std::cout << "  positive results: " << pos_results << std::endl;
        std::cout << "  negative results: " << neg_results << std::endl;
        std::cout << "  degenerate failures: " << degen_failures << std::endl;
        std::cout << std::endl;
    }
    std::cout << "  iteration count histogram:\n";
    std::cout << "    ";
    std::cout << "(max " << max_iters << ") ";
    for (index_t i = 0; i < iter_hist.size(); ++i) {
        std::cout << iter_hist[i] << " ";
    }
    std::cout << "\n";
    EXPECT_EQ(failures, 0);
}


constexpr index_t N_TESTS = 1'000'000;

TEST(TEST_MODULE_NAME, test_2f_gjk_sat_agree) {
    test_box_gjk_matches_sat<float, 2>(&rng, N_TESTS);
}

TEST(TEST_MODULE_NAME, test_2d_gjk_sat_agree) {
    test_box_gjk_matches_sat<double, 2>(&rng, N_TESTS);
}

TEST(TEST_MODULE_NAME, test_3f_gjk_sat_agree) {
    test_box_gjk_matches_sat<float, 3>(&rng, N_TESTS);
}

TEST(TEST_MODULE_NAME, test_3d_gjk_sat_agree) {
    test_box_gjk_matches_sat<double, 3>(&rng, N_TESTS);
}
