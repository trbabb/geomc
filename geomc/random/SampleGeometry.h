#pragma once

#include <geomc/random/DenseDistribution.h>
#include <geomc/shape/ShapeTypes.h>
#include <random>

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

/**
 * @defgroup random Random
 * @brief Sampling from random distributions.
 *
 * Classes conform to the patterns of the C++ `<random>` library and
 * are intercompatible with it.
 */

/// @addtogroup random
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
struct SampleShape {
    Shape shape;
    
    using shape_type  = Shape;
    using param_type  = Shape;
    using result_type = typename Shape::point_t;
    
    ShapeDistribution();
    ShapeDistribution(const Shape& s);
    
    void  reset();
    Shape param() const;
    void  param(const Shape& s);
    
    result_type min() requires BoundedObject<Shape>;
    result_type max() requires BoundedObject<Shape>;
    bool operator==(const ShapeDistribution& other) const;
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
        DenseUniformDistribution<T> u(0, 1);
        for (index_t i = 0; i < N + 1; ++i) {
            T xi = u(rng);
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
 */
template <typename T, index_t N>
struct SampleShape<Rect<T,N>> : public detail::ShapeDistribution<Rect<T,N>> {
    using detail::ShapeDistribution<Rect<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        using ptype = PointType<T,N>;
        if constexpr (N == 1) {
            DenseUniformDistribution<T> u {shape.lo, shape.hi};
            return u(rng);
        } else {
            DenseUniformDistribution<T> u(0, 1);
            point_t p;
            for (index_t i = 0; i < N; ++i) {
                T s = u(rng);
                ptype::iterator(p)[i] = s * shape.lo[i] + (1 - s) * shape.hi[i];
            }
            return p;
        }
    }
    
    bool operator==(const SampleShape& other) const = default;
};


/**
 * @brief Sample a point from the interior of a sphere.
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
            DenseUniformDistribution<T> unif {-1,1};
            point_t out;
            do {
                // generate a point in the signed unit box
                for (index_t i = 0; i < N; ++i) {
                    ptype::iterator(out)[i] = unif(rng);
                }
                // if the point is outside the unit sphere, reject it and try again
            } while (ptype::mag2(out) > 1);
            return shape.r * out + shape.center;
        } else {
            // draw a multivariate gaussian, then project it onto the sphere
            DenseUniformDistribution<T> u_01(0,1);
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
            DenseUniformDistribution<T> unif {-1,1};
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
        DenseUniformDistribution<T> u_dist(r0, r1);
        point_t p;
        if constexpr (N <= 3) {
            // rejection sampling
            DenseUniformDistribution<T> unif {-1,1};
            do {
                // generate a point in the signed unit box
                for (index_t i = 0; i < N; ++i) {
                    ptype::iterator(p)[i] = unif(rng);
                }
                // if the point is outside the unit sphere, reject it and try again
            } while (ptype::mag2(p) > 1);
        } else {
            // draw a multivariate gaussian, then project it onto the sphere
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
 */
template <typename T, index_t N>
struct SampleShape<Cylinder<T,N>> : public detail::ShapeDistribution<Cylinder<T,N>> {
    using detail::ShapeDistribution<Cylinder<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;

private:
    SampleShape<Sphere<T,N-1>> _sphere;
public:
    
    using detail::ShapeDistribution<Cylinder<T,N>>::ShapeDistribution;
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        constexpr point_t _Z = {VecType<T,N-1>((T)0), (T)1};
        DenseUniformDistribution<T> unif {0,1};
        // pick a point on the cap
        Vec<T,N> base {_sphere(rng) * shape.radius, 0};
        point_t axis = shape.p1 - shape.p0;
        // reorient the cap distribution to be perpendicular to the axis
        base = base.align(_Z, axis.unit());
        // pick a distance along the axis direction and follow it away from the base
        return base + (unif(rng) * axis) + shape.p0;
    }
    
    void reset() {
        _sphere.reset();
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/**
 * @brief Sample a point from the surface of a spherical cap.
 */
template <typename T, index_t N>
requires (N == 2 or N == 3)
struct SampleShape<SphericalCap<T,N>> : public detail::ShapeDistribution<SphericalCap<T,N>> {
    using detail::ShapeDistribution<SphericalCap<T,N>>::shape;
    using typename Dimensional<T,N>::point_t;
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        DenseUniformDistribution<T> uniform {0, 1};
        constexpr T pi = std::numbers::pi_v<T>;
        T half_angle = std::clamp<T>(shape.half_angle_radians, 0, pi);
        if constexpr (N == 2) {
            T u = uniform(rng);
            // pick theta in [-half_angle, +half_angle]
            T theta = u * half_angle - (1 - u) * half_angle;
            return {
                std::sin(theta),
                std::cos(theta)
            };
        } else if constexpr (N == 3) {
            T u = uniform(rng);
            T v = uniform(rng);
            // we pick a height with uniform probability
            // (because a slice through the sphere with thickness t at height h
            // has unchanging area regardless of h). 
            T min_h = std::cos(half_angle);
            T h = u * min_h + (1 - u);
            T r = std::sqrt(1 - h * h);
            T theta = v * 2 * pi;
            return {
                r * std::cos(theta),
                r * std::sin(theta),
                h
            };
        }
    }
    
};


/**
 * @brief Sample a point from the interior of a transformed shape.
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
        DenseUniformDistribution<T> unif {shape.height.lo, shape.height.hi};
        return { _inner(rng), unif(rng) };
    }
    
    bool operator==(const SampleShape& other) const {
        return shape == other.shape;
    }
};


/// @brief Sample a point from the boundary of a Rect.
template <typename T, index_t N>
struct SampleShape<Hollow<Rect<T,N>>> : public detail::ShapeDistribution<Hollow<Rect<T,N>>> {
    using detail::ShapeDistribution<Hollow<Rect<T,N>>>::shape;
    using typename Dimensional<T,N>::point_t;

private:
    point_t _face_areas;
    T       _face_area_sum;
public:
    
    SampleShape<Hollow<Rect<T,N>>>()        { param(shape); }
    SampleShape(const Hollow<Rect<T,N>>& s) { param(s); }
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        using ptype = PointType<T,N>;
        DenseUniformDistribution<T> u(0, 1);
        // choose a point inside the box
        point_t p;
        for (index_t i = 0; i < N; ++i) {
            T s = u(rng);
            ptype::iterator(p)[i] = s * shape.lo[i] + (1 - s) * shape.hi[i];
        }
        // choose a face in proportion to its area
        T sum = 0;
        index_t side = 0;
        T s = u(rng) * _face_area_sum;
        while (sum < s and side < N) {
            sum += _face_areas[side];
            ++side;
        }
        // project the point onto the chosen face
        T midpt = shape.lo[side] / 2 + shape.hi[side] / 2; // avoid overflow
        p[side] = p[side] > midpt ? shape.hi[side] : shape.lo[side];
        return p;
    }
    
    void param(const Hollow<Rect<T,N>>& s) {
        shape = s;
        _face_areas    = shape.face_areas();
        _face_area_sum = _face_areas.sum();
    }
    
    bool operator==(const SampleShape& other) const = default;
};


/// @brief Sample a point from the boundary of an Extruded shape.
template <typename Shape>
struct SampleShape<Hollow<Extruded<Shape>>> : public detail::ShapeDistribution<Hollow<Extruded<Shape>>> {
    using detail::ShapeDistribution<Hollow<Extruded<Shape>>>::shape;
    using typename detail::ShapeDistribution<Hollow<Extruded<Shape>>>::shape_type;
    using typename Dimensional<typename shape_type::elem_t, shape_type::N>::point_t;
    using T = typename shape_type::elem_t;
    using Dimensional<T,shape_type::N>::N;

private:
    SampleShape<Shape>         _face_shape_sampler;
    SampleShape<Hollow<Shape>> _wall_sampler;
    T _cap_area;
    T _wall_area;
public:
    
    SampleShape<Hollow<Extruded<Shape>>>()        { param(shape); }
    SampleShape(const Hollow<Extruded<Shape>>& s) { param(s); }
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        DenseUniformDistribution<T> u(0, 1);
        const Extruded<Shape>& extrusion = shape.shape;
        T s = u(rng) * (_cap_area + _wall_area);
        T h = u(rng);
        if (s < _cap_area) {
            // choose a point on the cap
            return {
                _face_shape_sampler(rng),
                h > 0.5 ? extrusion.height.hi : extrusion.height.lo
            };
        } else {
            // choose a point on the wall
            return {
                _wall_sampler(rng),
                extrusion.height.remap(h)
            };
        }
    }
    
    void param(const Hollow<Extruded<Shape>>& s) {
        shape = s;
        const Shape& face_shape = s.shape.base;
        _face_shape_sampler.param(face_shape);
        _wall_sampler.param(face_shape);
        _cap_area  = face_shape.measure_interior();
        _wall_area = face_shape.measure_boundary() * s.shape.height.dimensions();
    }
    
    bool operator==(const SampleShape& other) const = default;
};

/// @brief Sample a point from the surface of an extruded shape with no endcaps.
template <typename Shape>
struct SampleShape<Hollow<Extruded<Hollow<Shape>>>> :
        public detail::ShapeDistribution<Hollow<Extruded<Hollow<Shape>>>>
{
    using detail::ShapeDistribution<Hollow<Extruded<Hollow<Shape>>>>::shape;
    using typename detail::ShapeDistribution<Hollow<Extruded<Hollow<Shape>>>>::shape_type;
    using typename Dimensional<typename shape_type::elem_t, shape_type::N>::point_t;
    using T = typename shape_type::elem_t;
    using Dimensional<T,shape_type::N>::N;
    using face_shape_t = Hollow<Shape>;

private:
    SampleShape<face_shape_t> _face_shape_sampler;
public:
    
    SampleShape<shape_type>()        { param(shape); }
    SampleShape(const shape_type& s) { param(s); }
    
    template <typename Generator>
    point_t operator()(Generator& rng) {
        DenseUniformDistribution<T> u(0, 1);
        const Extruded<Shape>& extrusion = shape.shape;
        T h = u(rng);
        return {
            face_shape_sampler(rng),
            extrusion.height.remap(h)
        };
    }
    
    void param(const shape_type& s) {
        shape = s;
        const face_shape_t& face_shape = s.shape.base;
        _face_shape_sampler.param(face_shape);
    }
    
    bool operator==(const SampleShape& other) const = default;
};

/// @}

} // namespace geom
