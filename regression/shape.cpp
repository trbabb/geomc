#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Shape

// #include <iostream>

#include <random>
#include <boost/test/unit_test.hpp>
#include <geomc/function/Utils.h>
#include <geomc/shape/Oriented.h>
#include <geomc/shape/Cylinder.h>
#include <geomc/shape/Simplex.h>
#include <geomc/shape/Sphere.h>
#include <geomc/shape/Extrusion.h>
#include <geomc/shape/Frustum.h>


using namespace geom;
using namespace std;


typedef std::mt19937_64 rng_t;

// todo: test all the compound shapes.

// todo: tests:
//   - test convex_support against all the cardinal directions
//   - move a small distance away from a support point, both toward and away
//     from the object, and check for expected result of shape.contains()
//   - test that op(xf * shape, p) == op(shape, p / xf) for all {xf, p, shape, op}


rng_t rng(18374691138699945602ULL);


/****************************
 * random number generation *
 ****************************/


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
        Vec<T,N> p;
        for (index_t i = 0; i < N; ++i) {
            p[i] = shape.lo[i] + u(*rng) * (shape.hi[i] - shape.lo[i]);
        }
        return p;
    }
};


template <typename T, index_t N>
struct ShapeSampler<Sphere<T,N>> {
    Sphere<T,N> shape;
    ShapeSampler(const Sphere<T,N>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        std::uniform_real_distribution<T> u(0,1);
        Vec<T,N> p = rnd<T,N>(rng).unit();
        p = p * shape.r * std::pow(u(*rng), 1/(T)N) + shape.center;
        return p;
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
        Vec<T,N-1> px = ShapeSampler<Sphere<T,N-1>>(Sphere<T,N-1>(shape.radius))(rng);
        Vec<T,N> p = mix(u(*rng), shape.p0, shape.p1);
        for (index_t i = 0; i < N - 1; ++i) {
            p += bases[i + 1] * px[i];
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
        std::uniform_real_distribution<T> u(0,1);
        Vec<T,N-1> p = ShapeSampler<Shape>(shape.base)(rng);
        T h = (shape.height.hi - shape.height.lo) * u(*rng) + shape.height.lo;
        return Vec<T,N>(p, h);
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
        Vec<T,N> tx = 10* rnd<T,N>(rng);
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
        Vec<T,N>    n = rnd<T,N>(rng);
        Plane<T,N> pl = Plane<T,N>(n, s.convex_support(n));
        BOOST_CHECK(pl.contains(p));
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


template <template <typename> class Outer, typename T>
void explore_compound_shape(rng_t* rng, index_t shapes) {
    explore_shape<Outer<Rect<T, 2>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 3>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 4>>>(rng, shapes);
    explore_shape<Outer<Rect<T, 5>>>(rng, shapes);
    
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

BOOST_AUTO_TEST_SUITE_END()
