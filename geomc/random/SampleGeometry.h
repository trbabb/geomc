#pragma once

#include <random>
#include <geomc/shape/ShapeTypes.h>

namespace geom {
namespace detail {

template <typename Shape>
struct ShapeDistribution : public Dimensional<typename Shape::elem_t, Shape::N> {
    Shape shape;
    
    using shape_type  = Shape;
    using param_type  = Shape;
    using result_type = typename Shape::point_t;
    
    ShapeDistribution() = default;
    ShapeDistribution(const Shape& s):shape(s) {};
    
    void  reset() {}
    Shape param() const         { return shape; }
    void  param(const Shape& s) { shape = s; }
    
    result_type min() requires BoundedObject<Shape> {
        return shape.bounds().lo;
    }
    
    result_type max() requires BoundedObject<Shape> {
        return shape.bounds().hi;
    }
    
    bool operator==(const ShapeDistribution& other) const = default;
};

} // namespace detail

/// @ingroup random
/// @{

/**
 * @brief Sample a point from a shape.
 *
 * Points inside the shape have uniform probability of being sampled.
 * Use `SampleShape<Hollow<Shape>>` to sample points only on the boundary of a shape.
 *
 * Specializations only exist for certain shapes.
 */
#ifdef PARSING_DOXYGEN
template <typename Shape>
struct SampleShape : public detail::ShapeDistribution<Shape> {
    using typename Dimensional<typename Shape::elem_t, Shape::N>::point_t;
    
    /// Draw a point from the interior of the shape.
    template <typename Generator>
    point_t operator()(Generator& rng);
    
    bool operator==(const SampleShape& other) const;
};
#else
template <typename Shape>
struct SampleShape {};
#endif


/**
 * @brief Sample a point from the interior of a simplex.
 * @ingroup random
 */
template <typename T, index_t N>
struct SampleShape<Simplex<T,N>> : public detail::ShapeDistribution<Simplex<T,N>> {
    using detail::ShapeDistribution<Simplex<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;
    
private:
    std::exponential_distribution<T> e = 1;
public:
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        // algorithm from: https://projecteuclid.org/download/pdf_1/euclid.pjm/1102911301
        // we pick barycentric coordinates and then normalize so they sum to 1.
        T s[N + 1];
        T sum = 0;
        for (index_t i = 0; i < N + 1; ++i) {
            T xi = e(rng);
            s[i] = xi;
            sum += xi;
        }
        point_t p;
        for (index_t i = 0; i < N + 1; ++i) {
            T xi = s[i] / sum;
            p += xi * shape.pts[i];
        }
        return p;
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief @brief Sample a point from the interior of a rect.
 * @ingroup random
 */
template <typename T, index_t N>
struct SampleShape<Rect<T,N>> : public detail::ShapeDistribution<Rect<T,N>> {
    using detail::ShapeDistribution<Rect<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        using ptype = PointType<T,N>;
        std::uniform_real_distribution<T> u(0, 1);
        point_t p;
        for (index_t i = 0; i < N; ++i) {
            ptype::iterator(p)[i] = shape.lo[i] + u(rng) * (shape.hi[i] - shape.lo[i]);
        }
        return p;
    }
    
    bool operator==(const SampleShape& other) const = default;
};


/**
 * @brief Sample a point from the interior of a sphere.
 * @ingroup random
 */
template <typename T, index_t N>
struct SampleShape<Sphere<T,N>> : public detail::ShapeDistribution<Sphere<T,N>> {
    using detail::ShapeDistribution<Sphere<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;

private:
    using ptype = PointType<T,N>;
    using gauss_t = std::conditional_t<(N <= 3), int[0], std::normal_distribution<T>>;
    gauss_t _gauss;
public:
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        if constexpr (N <= 3) {
            // rejection sampling
            std::uniform_real_distribution<T> unif {-1,1};
            point_t out;
            do {
                // generate a point in the signed unit box
                for (index_t i = 0; i < N; ++i) {
                    ptype::iterator(out)[i] = unif(rng);
                }
                // if the point is outside the unit sphere, reject it and try again
            } while (ptype::mag2(out) > 1);
            return shape.radius * out + shape.center;
        } else {
            // draw a multivariate gaussian, then project it onto the sphere
            std::uniform_real_distribution<T> u_01(0,1);
            Vec<T,N> p;
            for (index_t i = 0; i < N; ++i) {
                p[i] = _gauss(rng);
            }
            return p.unit() * shape.r * std::pow(u_01(rng), 1 / (T)N) + shape.center;
        }
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from the surface of a sphere.
 * @ingroup random
 */
template <typename T, index_t N>
struct SampleShape<Hollow<Sphere<T,N>>> : public detail::ShapeDistribution<Hollow<Sphere<T,N>>> {
    using detail::ShapeDistribution<Hollow<Sphere<T,N>>>::shape;
    using typename Dimensional<T,N>::point_t;

private:
    using ptype = PointType<T,N>;
    using gauss_t = std::conditional_t<(N <= 3), int[0], std::normal_distribution<T>>;
    gauss_t _gauss;
public:
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        if constexpr (N <= 3) {
            // rejection sampling
            std::uniform_real_distribution<T> unif {-1,1};
            point_t out;
            do {
                // generate a point in the signed unit box
                for (index_t i = 0; i < N; ++i) {
                    ptype::iterator(out)[i] = unif(rng);
                }
                // if the point is outside the unit sphere, reject it and try again
            } while (ptype::mag2(out) > 1);
            return shape.radius * out.unit() + shape.center;
        } else {
            // draw a multivariate gaussian, then project it onto the sphere
            Vec<T,N> p;
            for (index_t i = 0; i < N; ++i) {
                p[i] = _gauss(rng);
            }
            return p.unit() * shape.r + shape.center;
        }
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from within a spherical shell.
 * @ingroup random
 */
template <typename T, index_t N>
struct SampleShape<SphericalShell<T,N>> : public detail::ShapeDistribution<SphericalShell<T,N>> {
    using detail::ShapeDistribution<SphericalShell<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;

private:
    using ptype = PointType<T,N>;
    using gauss_t = std::conditional_t<(N <= 3), int[0], std::normal_distribution<T>>;
    gauss_t _gauss;
public:
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        // recall that SphericalShell is Dilated<Hollow<Sphere>>
        //                         dilated.hollow.base
        const Sphere<T,N>& sphere = shape.shape.shape;
        T r0 = std::max<T>(sphere.radius - shape.dilation, 0);
        T r1 = sphere.radius + shape.dilation;
        std::uniform_real_distribution<T> u_dist(r0, r1);
        point_t p;
        if constexpr (N <= 3) {
            // rejection sampling
            std::uniform_real_distribution<T> unif {-1,1};
            do {
                // generate a point in the signed unit box
                for (index_t i = 0; i < N; ++i) {
                    ptype::iterator(p)[i] = unif(rng);
                }
                // if the point is outside the unit sphere, reject it and try again
            } while (ptype::mag2(p) > 1);
        } else {
            // draw a multivariate gaussian, then project it onto the sphere
            Vec<T,N> p;
            for (index_t i = 0; i < N; ++i) {
                p[i] = _gauss(rng);
            }
        }
        p = p.unit();
        T r = std::pow(u_dist(rng), 1 / (T)N);
        return p * r + sphere.center;
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from the interior of a cylinder.
 * @ingroup random
 */
template <typename T, index_t N>
struct SampleShape<Cylinder<T,N>> : public detail::ShapeDistribution<Cylinder<T,N>> {
    using detail::ShapeDistribution<Cylinder<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;

private:
    SampleShape<Sphere<T,N-1>> _sphere;
public:
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        constexpr point_t _Z = {{(T) 0}, (T) 1};
        std::uniform_real_distribution<T> unif {0,1};
        // pick a point on the cap
        auto base = _sphere(rng) * shape.radius;
        point_t axis = shape.p1 - shape.p0;
        // reorient the cap distribution to be perpendicular to the axis
        base = base.align(_Z, axis.unit());
        // pick a distance along the axis direction and follow it away from the base
        return base + unif(rng) * axis + shape.p0;
    }
    
    void reset() {
        _sphere.reset();
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from the interior of a transformed shape.
 * @ingroup random
 */
template <typename Shape>
struct SampleShape<Transformed<Shape>> : public detail::ShapeDistribution<Transformed<Shape>> {
    using detail::ShapeDistribution<Transformed<Shape>>::shape;
    
private:
    SampleShape<Shape> _inner;
public:
    
    SampleShape() = default;
    SampleShape(const SampleShape<Transformed<Shape>>& other) = default;
    SampleShape(SampleShape<Transformed<Shape>>&& other) = default;
    SampleShape(const Transformed<Shape>& s):
        detail::ShapeDistribution<Transformed<Shape>>(s),
        _inner(s.shape) {}
    
    void reset() {
        _inner.reset();
    }
    
    void param(const Transformed<Shape>& s) {
        shape = s;
        _inner.param(s.shape);
    }
    
    template <typename Generator>
    typename Shape::point_t operator()(Generator& rng) {
        return shape.xf * _inner(rng);
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from the interior of a shape transformed by a Similarity.
 * @ingroup random
 */
template <typename Shape>
struct SampleShape<Similar<Shape>> : public detail::ShapeDistribution<Similar<Shape>> {
    using detail::ShapeDistribution<Similar<Shape>>::shape;
    
private:
    SampleShape<Shape> _inner;
public:
    
    SampleShape() = default;
    SampleShape(const SampleShape<Similar<Shape>>& other) = default;
    SampleShape(SampleShape<Similar<Shape>>&& other) = default;
    SampleShape(const Similar<Shape>& s):
        detail::ShapeDistribution<Similar<Shape>>(s),
        _inner(s.shape) {}
    
    void reset() {
        _inner.reset();
    }
    
    void param(const Similar<Shape>& s) {
        shape = s;
        _inner.param(s.shape);
    }
    
    template <typename Generator>
    typename Shape::point_t operator()(Generator& rng) {
        return shape.xf * _inner(rng);
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from the interior of an extruded shape.
 * @ingroup random
 */
template <typename Shape>
struct SampleShape<Extruded<Shape>> : public detail::ShapeDistribution<Extruded<Shape>> {
    using detail::ShapeDistribution<Extruded<Shape>>::shape;
    
private:
    SampleShape<Shape> _inner;
public:
    
    SampleShape() = default;
    SampleShape(const SampleShape<Extruded<Shape>>& other) = default;
    SampleShape(SampleShape<Extruded<Shape>>&& other) = default;
    SampleShape(const Extruded<Shape>& s):
        detail::ShapeDistribution<Extruded<Shape>>(s),
        _inner(s.shape) {}
    
    void reset() {
        _inner.reset();
    }
    
    void param(const Extruded<Shape>& s) {
        shape = s;
        _inner.param(s.shape);
    }
    
    template <typename Generator>
    typename Shape::point_t operator()(Generator& rng) {
        using T = typename Shape::elem_t;
        std::uniform_real_distribution<T> unif {shape.height.lo, shape.height.hi};
        return { _inner(rng), unif(rng) };
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};

/// @}

} // namespace geom
