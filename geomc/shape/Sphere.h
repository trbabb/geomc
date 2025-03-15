#pragma once

#include <numbers>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Similarity.h>
#include <geomc/shape/Shape.h>
#include <geomc/function/Utils.h>

// todo: r should be "radius"

namespace geom {
   
/**
 * @brief Formula for the volume of a `d` dimensional ball with radius `r`.
 */
template <typename T>
constexpr T measure_ball_interior(index_t d, T r) {
    if (d == 0) return 1;
    if (d == 1) return 2 * r;
    return 2 * std::numbers::pi_v<T> * r * r * measure_ball_interior(d - 2, r) / (T) d;
}

/**
 * @brief Formula for the surface area of a `d` dimensional ball with radius `r`.
 */
template <typename T>
constexpr T measure_sphere_boundary(index_t d, T r) {
    constexpr T k = 2 * std::numbers::pi_v<T>;
    if (d == 0) return 1;
    if (d == 1) return 2;
    if (d == 2) return k * r;
    return k * r * r * measure_sphere_boundary(d - 2, r) / (T) (d - 2);
}

/** 
 * @ingroup shape
 * @brief An N-dimensional circle, sphere, or hypersphere with a filled interior.
 * 
 * `Circle<T>` is a template alias for `Sphere<T,2>`.
 */
template <typename T, index_t N>
class Sphere:
    public Convex          <T,N,Sphere<T,N>>,
    public RayIntersectable<T,N,Sphere<T,N>>,
    public Projectable     <T,N,Sphere<T,N>>
{
    typedef PointType<T,N> ptype;
public:
    using typename Dimensional<T,N>::point_t;
    /// Center of the sphere.
    point_t center;
    /// Radius of the sphere.
    T radius;
    
    /**
     * Construct a sphere at the origin with radius 1.
     */
    constexpr Sphere():
        center((T)0),
        radius(1) {}
    
    /**
     * Construct a sphere with center at the origin, having radius `r`.
     * @param r Radius of spehre.
     */
    constexpr Sphere(T r):
        center((T)0),
        radius(r) {}
    
    /**
     * Construct a sphere with center at the point `c`, having radius `r`.
     * @param c Center of sphere.
     * @param r Radius of spehre.
     */
    constexpr Sphere(const point_t& c, T r):
        center(c),
        radius(r) {}
    
    static constexpr bool admits_cusps() { return false; }
    
    bool operator==(const Sphere& other) const {
        return center == other.center && radius == other.radius;
    }
    
    Rect<T,N> bounds() const {
        point_t rvec(radius);
        return Rect<T,N>(center - rvec, center + rvec);
    }
    
    /// Shape-point intersection test.
    inline bool contains(point_t p) const {
        return ptype::mag2(center - p) <= radius * radius;
    }
    
    /**
     * Sphere-sphere intersection test.
     * @param s Another sphere.
     * @return `true` if `s` overlaps with this sphere's volume, false otherwise.
     */
    bool intersects(Sphere s) const {
        return ptype::mag2(center - s.center) <= radius * radius;
    }
    
    bool intersects(const Rect<T,N>& rect) const {
        return rect.dist2(radius) <= radius * radius;
    }
    
    point_t convex_support(point_t d) const {
        return center + ptype::unit(d) * radius;
    }
    
    /// Signed distance function.
    inline T sdf(point_t p) const {
        return ptype::mag(p - center) - radius;
    }
    
    /**
     * Return the point `p` orthogonally projected onto the surface of the shape.
     */
    inline point_t project(point_t p) const {
        return ptype::unit(p - center) * radius + center;
    }
    
    /// Outward-facing direction.
    inline point_t normal(point_t p) const {
        return ptype::unit(p - center);
    }
    
    /// Shape-ray intersection test.
    Rect<T,1> intersect(const Ray<T,N>& ray) const {
        T r2 = radius * radius;
        Vec<T,N> dir = ray.direction;
        Vec<T,N> x0 = center - ray.origin;
        // solve for s such  that ||s * ray - ctr|| == radius
        T a = dir.mag2();
        T b = -2 * dir.dot(x0);
        T c = x0.mag2() - r2;
        T roots[2];
        if (quadratic_solve(roots, a, b, c)) {
            return Rect<T,1>::from_corners(roots[0], roots[1]);
        } else {
            // empty interval
            return Rect<T,1>();
        }
    }
    
    /**
     * @brief Measure the interior (volume) of the shape.
     *
     * If the sphere is 2D (a disk), this is the area of the disk.
     * If the sphere is 3D (a ball), this is the volume of the ball.
     * In higher dimensions, this is the hypervolume.
     */
    T measure_interior() const {
        return measure_ball_interior(N, radius);
    }
    
    /**
     * @brief Measure the boundary of the shape.
     *
     * If the sphere is 2D (a circle), this is the circumference of the circle.
     * If the sphere is 3D (a sphere), this is the surface area of the sphere.
     * In higher dimensions, this is the volume or hypervolume of the boundary.
     */
    T measure_boundary() const {
        return measure_sphere_boundary(N, radius);
    }
    
}; /* Sphere<T,N> */

/// @addtogroup shape
/// @{

/// @brief Transform a sphere by a similarity transform.
/// @related Sphere
/// @related Similarity
template <typename T, index_t N>
Sphere<T,N> operator*(const Similarity<T,N>& xf, const Sphere<T,N>& s) {
    return Sphere<T,N>(xf * s.center, s.radius * xf.sx);
}

/// @brief Inverse-transform a sphere by a similarity transform.
/// @related Sphere
/// @related Similarity
template <typename T, index_t N>
Sphere<T,N> operator/(const Sphere<T,N>& s, const Similarity<T,N>& xf) {
    return Sphere<T,N>(xf / s.center, s.radius / xf.sx);
}

/// @brief Transform a sphere by an isometry.
/// @related Sphere
/// @related Isometry
template <typename T, index_t N>
Sphere<T,N> operator*(const Isometry<T,N>& xf, const Sphere<T,N>& s) {
    return Sphere<T,N>(xf * s.center, s.radius);
}

/// @brief Inverse-transform a sphere by an isometry.
/// @related Sphere
/// @related Isometry
template <typename T, index_t N>
Sphere<T,N> operator/(const Sphere<T,N>& s, const Isometry<T,N>& xf) {
    return Sphere<T,N>(xf / s.center, s.radius);
}

/// @} // addtogroup shape

template <typename T, index_t N, typename H>
struct Digest<Sphere<T,N>, H> {
    H operator()(const Sphere<T,N>& s) const {
        H nonce = geom::truncated_constant<H>(0x948904c693a7ddeb, 0xd121a1f8ce15ac6c);
        return geom::hash_many<H>(nonce, s.center, s.radius);
    }
};

#ifdef GEOMC_USE_STREAMS

template <typename T, index_t N>
std::ostream& operator<<(std::ostream& os, const Sphere<T,N>& s) {
    os << "Sphere(" << s.center << "," << s.radius << ")";
    return os;
}

#endif

} /* namespace geom */


template <typename T, index_t N>
struct std::hash<geom::Sphere<T,N>> {
    size_t operator()(const geom::Sphere<T,N> &s) const {
        return geom::hash<geom::Sphere<T,N>, size_t>(s);
    }
};
