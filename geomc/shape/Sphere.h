#pragma once

#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Similarity.h>
#include <geomc/shape/Shape.h>
#include <geomc/function/Utils.h>

namespace geom {

/** 
 * @ingroup shape
 * @brief An N-dimensional circle, sphere, or hypersphere.
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
    typedef typename ptype::point_t point_t;
    /// Center of the sphere.
    point_t center;
    /// Radius of the sphere.
    T r;

    /**
     * Construct a sphere at the origin with radius 1.
     */
    constexpr Sphere():
        center((T)0),
        r(1) {}
    
    /**
     * Construct a sphere with center at the origin, having radius `r`.
     * @param r Radius of spehre.
     */
    constexpr Sphere(T r):
        center((T)0),
        r(r) {}
    
    /**
     * Construct a sphere with center at the point `c`, having radius `r`.
     * @param c Center of sphere.
     * @param r Radius of spehre.
     */
    constexpr Sphere(const point_t& c, T r):
        center(c),
        r(r) {}
    
    bool operator==(const Sphere& other) const {
        return center == other.center && r == other.r;
    }
    
    Rect<T,N> bounds() const {
        point_t rvec(r);
        return Rect<T,N>(center - rvec, center + rvec);
    }
    
    /// Shape-point intersection test.
    inline bool contains(point_t p) const {
        return ptype::mag2(center - p) <= r * r;
    }
    
    /**
     * Sphere-sphere intersection test.
     * @param s Another sphere.
     * @return `true` if `s` overlaps with this sphere's volume, false otherwise.
     */
    bool intersects(Sphere s) const {
        return ptype::mag2(center - s.center) <= r * r;
    }
    
    bool intersects(const Rect<T,N>& rect) const {
        return rect.dist2(r) <= r * r;
    }
    
    point_t convex_support(point_t d) const {
        return center + ptype::unit(d) * r;
    }
    
    /// Signed distance function.
    inline T sdf(point_t p) const {
        return ptype::mag(p - center) - r;
    }
    
    /**
     * Return the point `p` orthogonally projected onto the surface of the shape.
     */
    inline point_t project(point_t p) const {
        return ptype::unit(p - center) * r + center;
    }
    
    /// Outward-facing direction.
    inline point_t normal(point_t p) const {
        return ptype::unit(p - center);
    }
    
    /// Shape-ray intersection test.
    Rect<T,1> intersect(const Ray<T,N>& ray) const {
        T r2 = r * r;
        Vec<T,N> dir = ray.direction;
        Vec<T,N> x0 = center - ray.origin;
        // solve for s such  that ||s * ray - ctr|| == radius
        T a = dir.mag2();
        T b = -2 * dir.dot(x0);
        T c = x0.mag2() - r2;
        T roots[2];
        if (quadratic_solve(roots, a, b, c)) {
            return Rect<T,1>::spanning_corners(roots[0], roots[1]);
        } else {
            // empty interval
            return Rect<T,1>();
        }
    }
    
}; /* Sphere<T,N> */

/// @addtogroup shape
/// @{

/// @brief Transform a sphere by a similarity transform.
/// @related Sphere
/// @related Similarity
template <typename T, index_t N>
Sphere<T,N> operator*(const Similarity<T,N>& xf, const Sphere<T,N>& s) {
    return Sphere<T,N>(xf * s.center, s.r * xf.sx);
}

/// @brief Inverse-transform a sphere by a similarity transform.
/// @related Sphere
/// @related Similarity
template <typename T, index_t N>
Sphere<T,N> operator/(const Sphere<T,N>& s, const Similarity<T,N>& xf) {
    return Sphere<T,N>(xf / s.center, s.r / xf.sx);
}

/// @brief Transform a sphere by an isometry.
/// @related Sphere
/// @related Isometry
template <typename T, index_t N>
Sphere<T,N> operator*(const Isometry<T,N>& xf, const Sphere<T,N>& s) {
    return Sphere<T,N>(xf * s.center, s.r);
}

/// @brief Inverse-transform a sphere by an isometry.
/// @related Sphere
/// @related Isometry
template <typename T, index_t N>
Sphere<T,N> operator/(const Sphere<T,N>& s, const Isometry<T,N>& xf) {
    return Sphere<T,N>(xf / s.center, s.r);
}

/// @} // addtogroup shape

template <typename T, index_t N, typename H>
struct Digest<Sphere<T,N>, H> {
    H operator()(const Sphere<T,N>& s) const {
        H nonce = geom::truncated_constant<H>(0x948904c693a7ddeb, 0xd121a1f8ce15ac6c);
        return geom::hash_many<H>(nonce, s.center, s.r);
    }
};

} /* namespace geom */


template <typename T, index_t N>
struct std::hash<geom::Sphere<T,N>> {
    size_t operator()(const geom::Sphere<T,N> &s) const {
        return geom::hash<geom::Sphere<T,N>, size_t>(s);
    }
};
