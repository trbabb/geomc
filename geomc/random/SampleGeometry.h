#pragma once

#include <random>
#include <geomc/shape/ShapeTypes.h>

namespace geom {
    
// todo: fix me when we implement Vec<T,1> decay


/// Implements the RandomDistribution concept
template <typename Shape>
struct ShapeDistribution {
    typedef Shape::point_t result_type;
    typedef Shape          param_type;
    Shape shape;
    
    ShapeDistribution(const Shape& s):shape(s) {}
    
    inline void  reset() {}
    inline Shape param() const         { return shape; }
    inline void  param(const Shape& s) { shape = s; }
};


template <typename Shape>
struct SampleShape {};


template <typename T, index_t N>
struct SampleShape<Simplex<T,N>> : public ShapeDistribution<Simplex<T,N>> {
    Simplex<T,N> shape;
    
    SampleShape(const Simplex<T,N>& s):ShapeDistribution(s) {}
    
    template <typename Generator>
    Vec<T,N> operator()(Generator& rng) {
        // algorithm from: https://projecteuclid.org/download/pdf_1/euclid.pjm/1102911301
        // we pick barycentric coordinates and then normalize so they sum to 1.
        std::exponential_distribution<T> e(1);
        T s[N + 1];
        T sum = 0;
        for (index_t i = 0; i < N + 1; ++i) {
            T xi = e(rng);
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
struct SampleShape<Rect<T,N>> : public ShapeDistribution<Rect<T,N>> {
    Rect<T,N> shape;
    SampleShape(const Rect<T,N>& s):ShapeDistribution(s) {}
    
    template <typename Generator>
    Vec<T,N> operator()(Generator& rng) {
        std::uniform_real_distribution<T> u(0, 1);
        Vec<T,N> p;
        for (index_t i = 0; i < N; ++i) {
            p[i] = shape.lo[i] + u(rng) * (shape.hi[i] - shape.lo[i]);
        }
        return p;
    }
};


template <typename T, index_t N>
struct SampleShape<Sphere<T,N>> : public ShapeDistribution<Sphere<T,N>> {
    Sphere<T,N> shape;
    SampleShape(const Sphere<T,N>& s):ShapeDistribution(s) {}
    
    template <typename Generator>
    Vec<T,N> operator()(Generator& rng) {
        std::uniform_real_distribution<T> u(0,1);
        Vec<T,N> p = rnd<T,N>(rng).unit();
        p = p * shape.r * std::pow(u(rng), 1/(T)N) + shape.center;
        return p;
    }
};


} // namespace geom
