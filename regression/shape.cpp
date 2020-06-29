#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Shape

// #include <iostream>

#include <random>
#include <boost/test/unit_test.hpp>
#include <geomc/function/Utils.h>
#include <geomc/shape/Oriented.h>
#include <geomc/shape/Cylinder.h>
#include <geomc/shape/Simplex.h>


using namespace geom;
using namespace std;


typedef std::mt19937_64 rng_t;

// todo: test all the shapes.

// todo: tests:
//   - test that random convex_support() also test true for `contains()`.
//   - test that the interpolation of two random `convex_support()` pass `contains()`.
//   - test that all such points are contained by bounds().
//   - test that op(xf * shape, p) == op(shape, p / xf) for all {xf, p, shape, op}

template <typename T, index_t N>
Vec<T,N> rnd(rng_t* rng) {
    std::normal_distribution<T> gauss(0, 1);
    Vec<T,N> v;
    for (index_t i = 0; i < N; ++i) {
        v[i] = gauss(*rng);
    }
    return v;
}

template <typename T>
inline T rnd(rng_t* rng) {
    std::normal_distribution<T> gauss(0, 1);
    return gauss(*rng);
}


template <typename Shape>
struct RandomShape {};

template <typename T, index_t N>
struct RandomShape<Sphere<T,N>> {
    static Sphere<T,N> rnd_shape(rng_t* rng) {
        return Sphere<T,N>(10 * rnd<T,N>(rng), 5 * std::abs(rnd<T>(rng)));
    }
};

template <typename T, index_t N>
struct RandomShape<Cylinder<T,N>> {
    static Cylinder<T,N> rnd_shape(rng_t* rng) {
        return Cylinder<T,N>(
            10 * rnd<T,N>(rng),
            10 * rnd<T,N>(rng),
            5 * std::abs(rnd<T>(rng)));
    }
};

template <typename T, index_t N>
struct RandomShape<Rect<T,N>> {
    static Rect<T,N> rnd_shape(rng_t* rng) {
        return Rect<T,N>::spanning_corners(
            10 * rnd<T,N>(rng),
            10 * rnd<T,N>(rng)
        );
    }
};

template <typename T, index_t N>
struct RandomShape<Simplex<T,N>> {
    static Simplex<T,N> rnd_shape(rng_t* rng) {
        Simplex<T,N> s;
        for (index_t i = 0; i <= N; ++i) {
            s += 10 * rnd<T,N>(rng);
        }
        return s;
    }
};

template <typename Shape>
struct RandomShape<Oriented<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr index_t N = Shape::N;
    
    static Oriented<Shape> rnd_shape(rng_t* rng) {
        SimpleMatrix<T,N,N> mx;
        for (index_t i = 0; i < N * N; ++i) {
            mx.begin()[i] = 3 * rnd<T>(rng);
        }
        AffineTransform<T,N> xf = translation(10 * rnd<T,N>(rng)) * transformation(mx);
        return Oriented<Shape>(RandomShape<Shape>::rnd_shape(rng), xf);
    }
};

template <typename Shape>
struct RandomShape<Frustum<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr index_t N = Shape::N;
    
    static Frustum<Shape> rnd_shape(rng_t* rng) {
        return Frustum<Shape>(
            RandomShape<Shape>::rnd_shape(rng),
            Rect<T,N>::spanning_corners(
                5 * rnd<T>(rng),
                5 * rnd<T>(rng)
            ));
    }
};

template <typename Shape>
struct RandomShape<Extrusion<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr index_t N = Shape::N;
    
    static Extrusion<Shape> rnd_shape(rng_t* rng) {
        return Extrusion<Shape>(
                RandomShape<Shape>::rnd_shape(rng),
                Rect<T,N>::spanning_corners(
                    5 * rnd<T>(rng),
                    5 * rnd<T>(rng)
                )
            );
    }
};


// xxx: this fails because interpolating verticies with themselves
// results in precision weirdness and failures.
template <typename Shape>
void exercise_shape(rng_t* rng, const Shape& s, index_t trials) {
    typedef typename Shape::elem_t T;
    constexpr index_t N = Shape::N;
    const T eps = 1e-5;
    std::uniform_real_distribution<T> interval(eps, 1 - eps);
    Vec<T,N> p0 = s.convex_support(rnd<T,N>(rng).unit());
    for (index_t i = 0; i < trials; ++i) {
        Vec<T,N> p1 = s.convex_support(rnd<T,N>(rng).unit());
        Vec<T,N>  p = mix(interval(*rng), p0, p1);
        // a random pt inside the shape should be contained by it:
        BOOST_CHECK(s.contains(p));
        // the convex support pts themselves should also be contained by the shape:
        BOOST_CHECK(s.contains(p1));
        p0 = p1;
    }
}


template <typename Shape>
void explore_shape(rng_t* rng, index_t shapes) {
    for (index_t i = 0; i < shapes; ++i) {
        Shape s = RandomShape<Shape>::rnd_shape(rng);
        exercise_shape<Shape>(rng, s, 10);
    }
}

BOOST_AUTO_TEST_SUITE(shape)


BOOST_AUTO_TEST_CASE(validate_rect) {
    rng_t rng(17581355241LL);
    
    // explore_shape<Rect<double, 2>>(&rng, 10);
    // explore_shape<Rect<double, 3>>(&rng, 10);
    // explore_shape<Rect<double, 4>>(&rng, 10);
}


BOOST_AUTO_TEST_CASE(create_oriented_cylinder) {
    // make a null-transformed oriented cylinder.
    // (the cylinder defaults to unit radius and length along X)
    auto ocyl = Oriented<Cylinder<double,3>>(Cylinder<double,3>());
    // confirm that the Oriented delegates containment checking
    BOOST_CHECK(ocyl.contains(Vec3d(0.5, 0, 0)));
    // confirm the Oriented delegages convex_support
    BOOST_CHECK_EQUAL(ocyl.convex_support(Vec3d(0.1, 1, 0)), Vec3d(1, 1, 0));
    // rotate the thing 180 degrees
    ocyl *= rotation(Vec3d(0, 0, 1), M_PI);
    // confirm that wrapper applies the xf:
    BOOST_CHECK(ocyl.contains(Vec3d(-0.5, 0, 0)));
}

BOOST_AUTO_TEST_CASE(orient_simple_shape) {
    auto xf = translation(Vec3d(-5, 0, 0));
    // confirm the operator works and its return type is correct:
    Oriented<Cylinder<double,3>> ocyl = xf * Cylinder<double,3>();
    // confirm the created wrapper applies the xf:
    BOOST_CHECK(ocyl.contains(Vec3d(-4.5, 0, 0)));
    // verify inheritance
    Convex<double,3>* s = &ocyl;
    s->convex_support(Vec3d(0.2,0.4,0.1));
}

BOOST_AUTO_TEST_SUITE_END()
