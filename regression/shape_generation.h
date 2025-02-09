#pragma once

#include <random>
#include <pcg_random.hpp>

#include <geomc/function/Utils.h>
#include <geomc/shape/Transformed.h>
#include <geomc/shape/Cylinder.h>
#include <geomc/shape/Simplex.h>
#include <geomc/shape/Sphere.h>
#include <geomc/shape/Plane.h>
#include <geomc/shape/Extruded.h>
#include <geomc/shape/Frustum.h>
#include <geomc/shape/SphericalCap.h>
#include <geomc/random/SampleGeometry.h>

using namespace geom;

// todo: better to generate shapes near the origin, and then translate them.
//   for isolated shapes, we might want big shapes far away.
//   for intersection tests, we want shapes in predictable places so we can make them touch

// pcg_extras::seed_seq_from<std::random_device> seed_source;

using rng_t = pcg64;

rng_t rng(6839552994771556349ULL);


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

template <typename T>
void rnd(rng_t* rng, AffineTransform<T,3>* xf) {
    Vec<T,4> z = rnd<T,4>(rng).unit();
    *xf = rotation(Quat<T>{z});
}

template <typename T>
void rnd(rng_t* rng, AffineTransform<T,2>* xf) {
    std::uniform_real_distribution<T> angle_range{-M_PI, M_PI};
    *xf = rotation(angle_range(*rng));
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
    
    // todo: use a householder reflection
    
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
struct ShapeSampler<Extruded<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr size_t N = Shape::N + 1;
    
    Extruded<Shape> shape;
    ShapeSampler(const Extruded<Shape>& s):shape(s) {}
    
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
        T s = u(*rng);
        T v = std::pow(s, 1./N); // ???
          v = (c_h.lo < 0) ? (1 - v) : v;
        T h = (c_h.hi - c_h.lo) * v + c_h.lo;
        return Vec<T,N>(h * p, h);
    }
};


template <typename Shape>
struct ShapeSampler<Transformed<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr size_t N = Shape::N;
    
    Transformed<Shape> shape;
    ShapeSampler(const Transformed<Shape>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        return shape.xf * ShapeSampler<Shape>(shape.shape)(rng);
    }
};

template <typename T, index_t N>
struct ShapeSampler<SphericalCap<T,N>> {
    SphericalCap<T,N> shape;
    ShapeSampler(const SphericalCap<T,N>& s):shape(s) {}
    
    Vec<T,N> operator()(rng_t* rng) {
        SampleShape<SphericalCap<T,N>> sampler(shape);
        return sampler(*rng);
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

template <typename T, index_t N>
struct RandomShape<SphericalCap<T,N>> {
    static SphericalCap<T,N> rnd_shape(rng_t* rng) {
        constexpr T pi = std::numbers::pi_v<T>;
        std::uniform_real_distribution<T> u(0, pi);
        return {u(*rng)};
    }
};

template <typename Shape>
struct RandomShape<Transformed<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr index_t N = Shape::N;
    
    static Transformed<Shape> rnd_shape(rng_t* rng) {
        SimpleMatrix<T,N,N> mx;
        for (index_t i = 0; i < N * N; ++i) {
            mx.begin()[i] = 3 * rnd<T>(rng);
        }
        AffineTransform<T,N> xf = translation(10 * rnd<T,N>(rng)) * transformation(mx);
        return Transformed<Shape>(RandomShape<Shape>::rnd_shape(rng), xf);
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
struct RandomShape<Extruded<Shape>> {
    typedef typename Shape::elem_t T;
    static constexpr index_t N = Shape::N;
    
    static Extruded<Shape> rnd_shape(rng_t* rng) {
        return Extruded<Shape>(
                RandomShape<Shape>::rnd_shape(rng),
                Rect<T,1>::spanning_corners(
                    5 * rnd<T>(rng),
                    5 * rnd<T>(rng)
                )
            );
    }
};
