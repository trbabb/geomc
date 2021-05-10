#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Shape

// #include <iostream>

#include <random>
#include <chrono>

#include <boost/test/unit_test.hpp>
#include <geomc/function/Utils.h>
#include <geomc/shape/Oriented.h>
#include <geomc/shape/Cylinder.h>
#include <geomc/shape/Simplex.h>
#include <geomc/shape/Sphere.h>
#include <geomc/shape/Extrusion.h>
#include <geomc/shape/Frustum.h>

#include "shapesample.h"


using namespace geom;
using namespace std;
using namespace std::chrono;

// todo: tests:
//   - move a small distance away from a support point, both toward and away
//     from the object, and check for expected result of shape.contains()
//   - test that op(xf * shape, p) == op(shape, p / xf) for all {xf, p, shape, op}

// todo: test Frustum<Rect<T,1>>

rng_t rng(18374691138699945602ULL);


/****************************
 * shape calisthenics       *
 ****************************/


template <typename Shape, typename T, index_t N>
void validate_plane(const Shape& s, const Vec<T,N>& p, const Vec<T,N>& n) {
    Plane<T,N> pl = Plane<T,N>(n, s.convex_support(n));
    BOOST_CHECK(pl.contains(p));
}


template <typename Shape>
bool validate_point(
        rng_t* rng, 
        const Shape& s, 
        const Vec<typename Shape::elem_t,Shape::N>& p) {
    typedef typename Shape::elem_t T;
    constexpr index_t N = Shape::N;
    if (s.contains(p)) {
        // if the pt is in the shape, it should definitely be in the bbox.
        BOOST_CHECK(s.bounds().contains(p));
        // pick a random support direction. the point should be inside
        // the bounds of the resultant plane.
        validate_plane(s, p, rnd<T,N>(rng));
        
        // test all the cardinal axes as support directions. a valid point
        // should come back, and the plane through it should still contain p.
        // (sometimes cardinal axes can be perfectly aligned with faces
        // and this degeneracy might be handled poorly).
        Vec<T,N> a;
        for (index_t i = 0; i < N; ++i) {
            a[i] = 1;
            validate_plane(s, p, a);
            a[i] = -1;
            validate_plane(s, p, a);
            a[i] =  0;
        }
        return true;
    }
    return false;
}


template <typename Shape>
void exercise_shape(rng_t* rng, const Shape& s, index_t trials) {
    typedef typename Shape::elem_t T;
    constexpr index_t N = Shape::N;
    auto sampler   = ShapeSampler<Shape>(s);
    Rect<T,N> bbox = s.bounds();
    // shrink the blob so that more of it falls inside the bbox:
    Vec<T,N>  dims = bbox.dimensions() / 4;
    Vec<T,N>     c = bbox.center();
    
    for (index_t i = 0; i < trials; ++i) {
        // validate a point near the shape's bbox
        validate_point<Shape>(rng, s, rnd<T,N>(rng) * dims + c);
        // validate a point definitely in the shape
        BOOST_CHECK(validate_point<Shape>(rng, s, sampler(rng)));
        // validate a point near the origin
        validate_point<Shape>(rng, s, rnd<T,N>(rng));
    }
}


template <typename Shape>
void explore_shape(rng_t* rng, index_t shapes) {
    for (index_t i = 0; i < shapes; ++i) {
        Shape s = RandomShape<Shape>::rnd_shape(rng);
        exercise_shape<Shape>(rng, s, 50);
    }
}


template <typename T, index_t N>
void explore_simplex(rng_t* rng, index_t shapes) {
    for (index_t i = 0; i < shapes; ++i) {
        Simplex<T,N> splx = RandomShape<Simplex<T,N>>::rnd_shape(rng);
        ShapeSampler<Simplex<T,N>> smp(splx);
        for (index_t j = 0; j < 100; ++j) {
            Vec<T,N> p = smp(rng);
            BOOST_CHECK(splx.contains(p));
        }
    }
}

template <typename T, index_t N>
void test_simplex_projection(rng_t* rng, const Simplex<T,N>& s, index_t trials) {
    auto bb   = s.bounds();
    auto dims = bb.dimensions();
    auto ctr  = bb.center();
    // no degenerate boxes pls
    for (index_t i = 0; i < N; ++i) {
        if (dims[i] == 0) dims[i] = 1;
    }
    for (index_t i = 0; i < trials; ++i) {
        Simplex<T,N> ss;
        Vec<T,N> p  = rnd<T,N>(rng) * dims + ctr;
        Vec<T,N> pp = s.project(p, &ss);
        if (ss.n == N + 1) {
            // `p` is inside a full-volume simplex.
            // no points should have been excluded; 
            // the "projected-to" simplex should be the same as the original
            BOOST_CHECK(ss == s);
            // if the point projected to the simplex volume, it should definitely 
            // be contained by the simplex
            BOOST_CHECK(ss.contains(p));
            // the "projected" point should be precisely the original point
            BOOST_CHECK_EQUAL(p, pp);
        } else {
            // the direction from the surface pt to the original point
            // should be orthogonal to the simplex face
            Vec<T,N> v = p - pp;
            for (index_t j = 1; j < ss.n; ++j) {
                Vec<T,N> b = ss.pts[j] - ss.pts[0];
                if (std::abs(b.dot(v)) > 1e-5) {
                    std::cout << ss.n << "/" << N << "\n";
                }
                BOOST_CHECK_SMALL(b.dot(v), 1e-5);
            }
            // the projected point should fall inside the face's bounds
            if (not ss.projection_contains(p)) {
                std::cout << "  " << ss.n << "/" << N << "\n";
            }
            BOOST_CHECK(ss.projection_contains(p));
        }
    }
}


template <typename T, index_t N>
void exercise_simplex_projection(rng_t* rng, index_t trials) {
    auto t0 = high_resolution_clock::now();
    for (index_t i = 0; i < trials; ++i) {
        Simplex<T,N> s = RandomShape<Simplex<T,N>>::rnd_shape(rng);
        std::uniform_int_distribution<> d(1,N+1);
        s.n = d(*rng); // make the simplex have a random number of verts
        test_simplex_projection(rng, s, 10);
    }
    auto t1 = high_resolution_clock::now();
    auto secs = duration_cast<duration<double>>(t1 - t0);
}


template <template <typename> class Outer, typename T>
void explore_compound_shape(rng_t* rng, index_t shapes) {
    explore_shape<Outer<Rect<T, 2>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 3>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 4>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 5>>>(rng, shapes);
    
    // explore_shape<Outer<Cylinder<T, 2>>>(rng, shapes); // [1]
    explore_shape<Outer<Cylinder<T, 3>>>(rng, shapes);
    explore_shape<Outer<Cylinder<T, 4>>>(rng, shapes);
    explore_shape<Outer<Cylinder<T, 5>>>(rng, shapes);
    explore_shape<Outer<Cylinder<T, 7>>>(rng, shapes);
    
    explore_shape<Outer<Sphere<T, 2>>>(rng, shapes);
    explore_shape<Outer<Sphere<T, 3>>>(rng, shapes);
    explore_shape<Outer<Sphere<T, 4>>>(rng, shapes);
    explore_shape<Outer<Sphere<T, 5>>>(rng, shapes);
    explore_shape<Outer<Sphere<T, 7>>>(rng, shapes);
    
    explore_shape<Outer<Simplex<T, 2>>>(rng, shapes);
    explore_shape<Outer<Simplex<T, 3>>>(rng, shapes);
    explore_shape<Outer<Simplex<T, 4>>>(rng, shapes);
    explore_shape<Outer<Simplex<T, 5>>>(rng, shapes);
    explore_shape<Outer<Simplex<T, 7>>>(rng, shapes);
    
    // [1] A valid construction, but not currently tested
    // because the sampling code for Cylinder needs to
    // draw from an N-1 sphere. When N=2, this is a line,
    // which degenerates to a different point_t, and Sphere
    // is not set up to work with N=1.
}


/****************************
 * test cases               *
 ****************************/


BOOST_AUTO_TEST_SUITE(shape)


BOOST_AUTO_TEST_CASE(validate_rect) {
    explore_shape<Rect<double, 2>>(&rng, 1000);
    explore_shape<Rect<double, 3>>(&rng, 1000);
    explore_shape<Rect<double, 4>>(&rng, 1000);
    explore_shape<Rect<double, 5>>(&rng, 1000);
}

BOOST_AUTO_TEST_CASE(validate_cylinder) {
    explore_shape<Cylinder<double, 3>>(&rng, 1000);
    explore_shape<Cylinder<double, 4>>(&rng, 1000);
    explore_shape<Cylinder<double, 5>>(&rng, 1000);
    explore_shape<Cylinder<double, 7>>(&rng, 1000);
}

BOOST_AUTO_TEST_CASE(validate_simplex) {
    explore_shape<Simplex<double, 2>>(&rng, 1000);
    explore_shape<Simplex<double, 3>>(&rng, 1000);
    explore_shape<Simplex<double, 4>>(&rng, 1000);
    explore_shape<Simplex<double, 5>>(&rng, 1000);
    explore_shape<Simplex<double, 7>>(&rng, 1000);
    // todo: also check that contains() and projection_contains(),
    //       all agree about pt containment.
}

BOOST_AUTO_TEST_CASE(validate_sphere) {
    explore_shape<Sphere<double, 2>>(&rng, 1000);
    explore_shape<Sphere<double, 3>>(&rng, 1000);
    explore_shape<Sphere<double, 4>>(&rng, 1000);
}


BOOST_AUTO_TEST_CASE(validate_extrusion) {
    explore_compound_shape<Extrusion, double>(&rng, 250);
}

BOOST_AUTO_TEST_CASE(validate_oriented) {
    explore_compound_shape<Oriented, double>(&rng, 250);
}

BOOST_AUTO_TEST_CASE(validate_frustum) {
    explore_compound_shape<Frustum, double>(&rng, 250);
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

BOOST_AUTO_TEST_CASE(simplex_projection) {
    const index_t iters = 10000;
    exercise_simplex_projection<double,2>(&rng, iters);
    exercise_simplex_projection<double,3>(&rng, iters);
    exercise_simplex_projection<double,4>(&rng, iters);
    exercise_simplex_projection<double,5>(&rng, iters);
    exercise_simplex_projection<double,7>(&rng, iters);
}

BOOST_AUTO_TEST_SUITE_END()
