/*
 * Sphere.h
 *
 *  Created on: May 10, 2009
 *      Author: Tim Babb
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include <geomc/linalg/Vec.h>
#include <geomc/shape/Shape.h>
#include <geomc/shape/shapedetail/Hit.h>
#include <geomc/function/Utils.h>

namespace geom {

/** 
 * @ingroup shape
 * @brief A N-dimensional circle, sphere, or hypersphere
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
    inline bool intersects(Sphere s) const {
        return ptype::mag2(center - s.center) <= r * r;
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

} /* namespace geom */
#endif /* SPHERE_H_ */
