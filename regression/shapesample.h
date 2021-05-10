#ifndef __SHAPESAMPLE_H_
#define __SHAPESAMPLE_H_

#include <random>

#include <geomc/shape/Oriented.h>
#include <geomc/shape/Cylinder.h>
#include <geomc/shape/Simplex.h>
#include <geomc/shape/Sphere.h>
#include <geomc/shape/Extrusion.h>
#include <geomc/shape/Frustum.h>

namespace geom {


typedef std::mt19937_64 rng_t;


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
        T v = std::pow(u(*rng), 1/(T)N);
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

};


#endif  // __SHAPESAMPLE_H_

