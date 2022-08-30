#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Shape

// #include <iostream>

#include <random>
#include <boost/test/unit_test.hpp>
#include <boost/core/demangle.hpp>

#include <geomc/function/Utils.h>
#include <geomc/shape/Oriented.h>
#include <geomc/shape/Cylinder.h>
#include <geomc/shape/Simplex.h>
#include <geomc/shape/Sphere.h>
#include <geomc/shape/Plane.h>
#include <geomc/shape/Extrusion.h>
#include <geomc/shape/Frustum.h>

// increase for good coverage.
// set to 1 for debugging.
#define N_TESTS 1

using namespace geom;
using namespace std;

typedef std::mt19937_64 rng_t;

// todo: tests:
//   - move a small distance away from a support point, both toward and away
//     from the object, and check for expected result of shape.contains()
//   - test that op(xf * shape, p) == op(shape, p / xf) for all {xf, p, shape, op}
//   - trace a ray and verify sdf(hit) << 1
//   - pick a segment that does not intersect the bbox of a shape
//     - verify no ray hit
//   - projection orthogonality:
//     - sample around the projected point. if any points satisfy:
//       - sampled pt inside the shape
//       - sampled pt closer to p than the projected point
//       - p outside the shape
//       ...then the projected point was not the closest; reject.

// todo: test Frustum<Rect<T,1>>
//   >> needs impl

rng_t rng(18374691138699945602ULL);


/****************************
 * random number generation *
 ****************************/


template <typename T, index_t N>
typename PointType<T,N>::point_t rnd(rng_t* rng) {
    std::normal_distribution<T> gauss(0, 1);
    typename PointType<T,N>::point_t v;
    for (index_t i = 0; i < N; ++i) {
        PointType<T,N>::iterator(v)[i] = gauss(*rng);
    }
    return v;
}

template <typename T>
inline T rnd(rng_t* rng) {
    std::normal_distribution<T> gauss(0, 1);
    return gauss(*rng);
}


/****************************
 * point sampling of shapes *
 ****************************/


template <typename Shape>
struct ShapeSampler {};


template <typename T, index_t N>
struct ShapeSampler<Simplex<T,N>> {
    Simplex<T,N> shape;
    
    ShapeSampler(const Simplex<T,N>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        // algorithm from: https://projecteuclid.org/download/pdf_1/euclid.pjm/1102911301
        // we pick barycentric coordinates and then normalize so they sum to 1.
        std::exponential_distribution<T> e(1);
        T s[N + 1];
        T sum = 0;
        for (index_t i = 0; i < N + 1; ++i) {
            T xi = e(*rng);
            s[i] = xi;
            sum += xi;
        }
        Vec<T,N> p;
        for (index_t i = 0; i < N + 1; ++i) {
            T xi = s[i] / sum;
            p += xi * shape.pts[i];
        }
        return p;
    }
};


template <typename T, index_t N>
struct ShapeSampler<Rect<T,N>> {
    Rect<T,N> shape;
    ShapeSampler(const Rect<T,N>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        std::uniform_real_distribution<T> u(0, 1);
        // pick a point in the unit cube:
        typename Rect<T,N>::point_t p;
        T* pp = PointType<T,N>::iterator(p);
        for (index_t i = 0; i < N; ++i) {
            pp[i] = u(*rng);
        }
        // map it to the shape:
        return shape.remap(p);
    }
};


template <typename T, index_t N>
struct ShapeSampler<Sphere<T,N>> {
    Sphere<T,N> shape;
    ShapeSampler(const Sphere<T,N>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        std::uniform_real_distribution<T> u(0,1);
        auto p = PointType<T,N>::unit(rnd<T,N>(rng));
        p = p * shape.r * std::pow(u(*rng), 1./N) + shape.center;
        return p;
    }
};

template <typename T>
struct ShapeSampler<Sphere<T,1>> {
    Sphere<T,1> shape;
    ShapeSampler(const Sphere<T,1>& s):shape(s) {}
    
    T operator()(rng_t* rng) {
        std::uniform_real_distribution<T> u(shape.center - shape.r, shape.center + shape.r);
        return u(*rng);
    }
};


template <typename T, index_t N>
struct ShapeSampler<Cylinder<T,N>> {
    Cylinder<T,N> shape;
    Vec<T,N> bases[N];
    
    ShapeSampler(const Cylinder<T,N>& s):shape(s) {
        bases[0] = s.p1 - s.p0;
        nullspace(bases, 1, bases + 1);
        orthonormalize(bases, N);
    }
    
    Vec<T,N> operator()(rng_t* rng) {
        std::uniform_real_distribution<T> u(0,1);
        auto px = ShapeSampler<Sphere<T,N-1>>(Sphere<T,N-1>(shape.radius))(rng);
        Vec<T,N> p = mix(u(*rng), shape.p0, shape.p1);
        for (index_t i = 0; i < N - 1; ++i) {
            T* w = PointType<T,N-1>::iterator(px);
            p += bases[i + 1] * w[i];
        }
        return p;
    }
};


template <typename Shape>
struct ShapeSampler<Extrusion<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr size_t N = Shape::N + 1;
    
    Extrusion<Shape> shape;
    ShapeSampler(const Extrusion<Shape>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        std::uniform_real_distribution<T> u(shape.height.lo, shape.height.hi);
        Vec<T,N-1> p = ShapeSampler<Shape>(shape.base)(rng);
        return Vec<T,N>(p, u(*rng));
    }
};


template <typename Shape>
struct ShapeSampler<Frustum<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr size_t N = Shape::N + 1;
    
    Frustum<Shape> shape;
    ShapeSampler(const Frustum<Shape>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        // i am not positive this logic is right, but at worst
        // we just have a skewed sampling of our frustum
        std::uniform_real_distribution<T> u(0,1);
        Vec<T,N-1> p = ShapeSampler<Shape>(shape.base)(rng);
        auto c_h = shape.clipped_height();
        T v = std::sqrt(u(*rng));
          v = (c_h.lo < 0) ? (1 - v) : v;
        T h = (c_h.hi - c_h.lo) * v + c_h.lo;
        return Vec<T,N>(h * p, h);
    }
};


template <typename Shape>
struct ShapeSampler<Oriented<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr size_t N = Shape::N;
    
    Oriented<Shape> shape;
    ShapeSampler(const Oriented<Shape>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        return shape.xf * ShapeSampler<Shape>(shape.shape)(rng);
    }
};


/****************************
 * random shape generation  *
 ****************************/


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
        // try not to be consistently centered on the origin:
        Vec<T,N> tx = 10 * rnd<T,N>(rng);
        return Cylinder<T,N>(
            5 * rnd<T,N>(rng) + tx,
            5 * rnd<T,N>(rng) + tx,
            2 * std::abs(rnd<T>(rng)));
    }
};

template <typename T, index_t N>
struct RandomShape<Rect<T,N>> {
    static Rect<T,N> rnd_shape(rng_t* rng) {
        auto tx = 10* rnd<T,N>(rng);
        return Rect<T,N>::spanning_corners(
            5 * rnd<T,N>(rng) + tx,
            5 * rnd<T,N>(rng) + tx
        );
    }
};

template <typename T, index_t N>
struct RandomShape<Simplex<T,N>> {
    static Simplex<T,N> rnd_shape(rng_t* rng) {
        Simplex<T,N> s;
        Vec<T,N> tx = 15 * rnd<T,N>(rng);
        for (index_t i = 0; i <= N; ++i) {
            s |= 5 * rnd<T,N>(rng) + tx;
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
            Rect<T,1>::spanning_corners(
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
                Rect<T,1>::spanning_corners(
                    5 * rnd<T>(rng),
                    5 * rnd<T>(rng)
                )
            );
    }
};


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
        const typename Shape::point_t& p) {
    typedef typename Shape::elem_t T;
    constexpr index_t N = Shape::N;
    if (s.contains(p)) {
        // if the pt is in the shape, it should definitely be in the bbox.
        BOOST_CHECK(s.bounds().contains(p));
        // pick a random support direction. the point should be inside
        // the bounds of the resultant plane.
        if constexpr (N > 1) {
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
        }
        return true;
    }
    return false;
}


template <typename Shape>
void validate_sdf(rng_t* rng, const Shape& s, index_t trials) {
    if constexpr (implements_shape_concept<Shape, SdfEvaluable>::value) {
        typedef typename Shape::elem_t T;
        constexpr index_t N = Shape::N;
        auto bb   = s.bounds();
        auto dims = bb.dimensions();
        auto ctr  = bb.center();
        for (index_t i = 0; i < trials; ++i) {
            auto p = rnd<T,N>(rng) * dims + ctr;
            BOOST_CHECK((s.sdf(p) <= 0) == s.contains(p));
        }
    }
}


template <typename Shape>
void validate_projection(rng_t* rng, const Shape& s, index_t trials) {
    if constexpr (implements_shape_concept<Shape, Projectable>::value) {
        typedef typename Shape::elem_t T;
        typedef typename Shape::point_t point_t;
        typedef PointType<T,Shape::N> ptype;
        constexpr index_t N = Shape::N;
        auto bb   = s.bounds();
        auto dims = bb.dimensions();
        auto ctr  = bb.center();
        for (index_t i = 0; i < trials; ++i) {
            auto p  = rnd<T,N>(rng) * dims + ctr;
            auto pp = s.project(p);
            T   sdf = s.sdf(pp);
            // projected point should be very near the surface
            BOOST_CHECK_SMALL(sdf, 1e-5);
            // projection should be idempotent
            // (in cases where sdf() has distinct implementation)
            T sz = ptype::mag(pp - s.project(pp));
            BOOST_CHECK_SMALL(sz, 1e-5);
            
            // the projection direction should be approximately normal to the surface
            point_t axis = p - pp;
            point_t n;
            T dL = ptype::mag(axis) * 1e-4;
            // compute the normal by taking the gradient of the sdf
            for (index_t j = 0; j < N; ++j) {
                auto dPdj = p;
                ptype::iterator(dPdj)[j] += dL;
                ptype::iterator(n)[j]     = (s.sdf(dPdj) - sdf) / dL;
            }
            // xxx failing for, like, every shape. probably a flaw with the test
            // BOOST_CHECK_SMALL(std::abs(n.unit().dot(axis.unit())) - 1, 1e-3);
            
            // a ray intersection should agree with project() about where the shape is
            if constexpr (implements_shape_concept<Shape, RayIntersectable>::value) {
                // cast a ray from the point along the projection direction.
                Ray<T,N> r = Ray<T,N>(p, -axis);
                Rect<T,1> interval = s.intersect(r);
                T t = interval.lo > 0 ? interval.lo : interval.hi;
                // the first positive hit should be at the projected point
                T raydist = ptype::mag(r.at_multiple(t) - pp);
                if (std::abs(raydist) > 1e-4) {
                    std::cout << p << " -> " << pp << "\n";
                }
                BOOST_CHECK_SMALL(raydist, 1e-4);
                // the hit should be at t=~1
                BOOST_CHECK_SMALL(t - 1, 1e-4);
            }
        }
    }
}


template <typename Shape>
void validate_ray(rng_t* rng, const Shape& s, index_t trials) {
    if constexpr (implements_shape_concept<Shape, RayIntersectable>::value) {
        typedef typename Shape::elem_t T;
        typedef typename Shape::point_t point_t;
        typedef PointType<T,Shape::N> ptype;
        constexpr index_t N = Shape::N;
        auto sampler = ShapeSampler<Shape>(s);
        for (index_t i = 0; i < trials; ++i) {
            // pick a point in the shape
            auto p = sampler(rng);
            // pick a random direction
            auto v = rnd<T,N>(rng);
            // construct a ray through the point
            Ray<T,N> ray{p,v};
            // intersect the ray with the shape
            auto interval = s.intersect(ray);
            // ray goes through a point in the shape. ray hit should reflect that:
            BOOST_CHECK(not interval.is_empty());
            // the ray origin is a point in the shape; s=0 should be in the overlap region
            BOOST_CHECK(interval.contains(0));
            // ray in the opposite direction should still intersect
            BOOST_CHECK(not s.intersect(Ray<T,N>{p, -v}).is_empty());
            // a ray displaced randomly along V should still intersect
            auto k = rnd<T>(rng);
            BOOST_CHECK(not s.intersect(Ray<T,N>{ray.at_multiple(k), v}).is_empty());
        }
    }
}


template <typename Shape>
void exercise_shape(rng_t* rng, const Shape& s, index_t trials) {
    typedef typename Shape::elem_t  T;
    typedef typename Shape::point_t point_t;
    constexpr index_t N = Shape::N;
    auto sampler   = ShapeSampler<Shape>(s);
    Rect<T,N> bbox = s.bounds();
    // shrink the blob so that more of it falls inside the bbox:
    point_t  dims = bbox.dimensions() / 4;
    point_t     c = bbox.center();
    
    for (index_t i = 0; i < trials; ++i) {
        // validate a point near the shape's bbox
        validate_point<Shape>(rng, s, rnd<T,N>(rng) * dims + c);
        // validate a point definitely in the shape
        BOOST_CHECK(validate_point<Shape>(rng, s, sampler(rng)));
        // validate a point near the origin
        validate_point<Shape>(rng, s, rnd<T,N>(rng));
    }
    
    validate_sdf(rng, s, trials);
    validate_projection(rng, s, trials);
    validate_ray(rng, s, trials);
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
                BOOST_CHECK_SMALL(b.dot(v), 1e-5);
            }
            // the projected point should fall inside the face's bounds
            BOOST_CHECK(ss.projection_contains(p));
        }
    }
}


template <typename T, index_t N>
void exercise_simplex_projection(rng_t* rng, index_t trials) {
    for (index_t i = 0; i < trials; ++i) {
        Simplex<T,N> s = RandomShape<Simplex<T,N>>::rnd_shape(rng);
        std::uniform_int_distribution<> d(1,N+1);
        s.n = d(*rng); // make the simplex have a random number of verts
        test_simplex_projection(rng, s, 10);
    }
}


template <template <typename> class Outer, typename T, bool test_degenerate=false>
void explore_compound_shape(rng_t* rng, index_t shapes) {
    explore_shape<Outer<Rect<T, 2>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 3>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 4>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 5>>>(rng, shapes);
    
    explore_shape<Outer<Cylinder<T, 2>>>(rng, shapes);
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
    
    if constexpr (test_degenerate) {
        explore_shape<Outer<Rect<T, 1>>>(rng, shapes);
        explore_shape<Outer<Sphere <T, 1>>>(rng, shapes);
        explore_shape<Outer<Simplex<T, 1>>>(rng, shapes);
    }
}


/****************************
 * test cases               *
 ****************************/


BOOST_AUTO_TEST_SUITE(shape)


BOOST_AUTO_TEST_CASE(validate_rect) {
    explore_shape<Rect<double, 2>>(&rng, N_TESTS);
    explore_shape<Rect<double, 3>>(&rng, N_TESTS);
    explore_shape<Rect<double, 4>>(&rng, N_TESTS);
    explore_shape<Rect<double, 5>>(&rng, N_TESTS);
}

BOOST_AUTO_TEST_CASE(validate_cylinder) {
    explore_shape<Cylinder<double, 2>>(&rng, N_TESTS);
    explore_shape<Cylinder<double, 3>>(&rng, N_TESTS);
    explore_shape<Cylinder<double, 4>>(&rng, N_TESTS);
    explore_shape<Cylinder<double, 5>>(&rng, N_TESTS);
    explore_shape<Cylinder<double, 7>>(&rng, N_TESTS);
}

BOOST_AUTO_TEST_CASE(validate_simplex) {
    explore_shape<Simplex<double, 2>>(&rng, N_TESTS);
    explore_shape<Simplex<double, 3>>(&rng, N_TESTS);
    explore_shape<Simplex<double, 4>>(&rng, N_TESTS);
    explore_shape<Simplex<double, 5>>(&rng, N_TESTS);
    explore_shape<Simplex<double, 7>>(&rng, N_TESTS);
    // todo: also check that contains() and projection_contains(),
    //       all agree about pt containment.
}

BOOST_AUTO_TEST_CASE(validate_sphere) {
    explore_shape<Sphere<double, 1>>(&rng, N_TESTS);
    explore_shape<Sphere<double, 2>>(&rng, N_TESTS);
    explore_shape<Sphere<double, 3>>(&rng, N_TESTS);
    explore_shape<Sphere<double, 4>>(&rng, N_TESTS);
}


BOOST_AUTO_TEST_CASE(validate_extrusion) {
    explore_compound_shape<Extrusion, double>(&rng, std::max(N_TESTS / 4, 1));
}

BOOST_AUTO_TEST_CASE(validate_oriented) {
    explore_compound_shape<Oriented, double>(&rng, std::max(N_TESTS / 4, 1));
}

BOOST_AUTO_TEST_CASE(validate_frustum) {
    // explore_compound_shape<Frustum, double>(&rng, std::max(N_TESTS / 4, 1));
    explore_shape<Frustum<Rect<double,2>>>(&rng, 1);
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
    ocyl.convex_support(Vec3d(0.2,0.4,0.1));
}

BOOST_AUTO_TEST_CASE(simplex_projection) {
    exercise_simplex_projection<double,2>(&rng, N_TESTS);
    exercise_simplex_projection<double,3>(&rng, N_TESTS);
    exercise_simplex_projection<double,4>(&rng, N_TESTS);
    exercise_simplex_projection<double,5>(&rng, N_TESTS);
    exercise_simplex_projection<double,7>(&rng, N_TESTS);
}

BOOST_AUTO_TEST_SUITE_END()
