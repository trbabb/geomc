#pragma once

#include <geomc/linalg/Vec.h>
#include <geomc/shape/Shape.h>
#include <geomc/function/Utils.h>

namespace geom {

/** 
 * @ingroup shape
 * @brief A N-dimensional capsule shape
 */
template <typename T, index_t N>
class Capsule:
    public Convex          <T,N,Capsule<T,N>>,
    public RayIntersectable<T,N,Capsule<T,N>>,
    public Projectable     <T,N,Capsule<T,N>>
{
    typedef PointType<T,N> ptype;
public:
    typedef typename ptype::point_t point_t;
    /// endpoints of the capsule axis.
    point_t p0;
    point_t p1;
    /// Radius of the capsule; i.e. the distance from the axis which
    /// is inside the capsule.
    T radius;

    /**
     * Construct a capsule centered on the origin, with radius 1, and axis
     * from `x = -1` to `x = 1`.
     */
    constexpr Capsule():radius(1) {
        p0[0] = -1;
        p1[0] =  1;
    }
    
    constexpr Capsule(point_t p0, point_t p1, T r):
        p0(p0),
        p1(p1),
        radius(r) {}
    
    bool operator==(const Capsule& other) const {
        return radius == other.radius and p0 == other.p0 and p1 == other.p1;
    }
    
    Rect<T,N> bounds() const {
        return Rect<T,N>::spanning_corners(p0, p1).dilated(radius);
    }
    
    /// Shape-point intersection test.
    bool contains(point_t p) const {
        return sdf(p) <= 0;
    }
    
    point_t nearest_axis_point(point_t p) const {
        point_t v    = p  - p0;
        point_t axis = p1 - p0;
        T frac = v.dot(axis) / axis.mag2();
        return p0 + axis * clamp(frac, 0, 1);
    }
    
    /**
     * Sphere-capsule intersection test.
     * @param s Another sphere.
     * @return `true` if `s` overlaps with this sphere's volume, false otherwise.
     */
    bool intersects(Sphere<T,N> s) const {
        return sdf(s.center) <= s.r;
    }
    
    point_t convex_support(point_t d) const {
        point_t axis = p1 - p0;
        point_t p = d.dot(axis) >= 0 ? p1 : p0;
        return p + d.with_length(radius);
    }
    
    /// Signed distance function.
    T sdf(point_t p) const {
        point_t b    = p  - p0;
        point_t axis = p1 - p0;
        T frac = clamp<T>(b.dot(axis) / axis.mag2(), 0, 1);
        return (b - frac * axis).mag() - radius;
    }
    
    /**
     * Return the point `p` orthogonally projected onto the surface of the shape.
     */
    point_t project(point_t p) const {
        point_t axis_pt = nearest_axis_point(p);
        point_t v = p - axis_pt;
        return axis_pt + v.with_length(radius);
    }
    
    /// Outward-facing direction.
    point_t normal(point_t p) const {
        return ptype::unit(p - nearest_axis_point(p));
    }
    
    /*
    xxx: todo: this. consider using https://www.shadertoy.com/view/Xt3SzX
    /// Shape-ray intersection test.
    Rect<T,1> intersect(const Ray<T,N>& ray) const {
        Vec<T,N> a = p1 - p0;
        Vec<T,N> b = ray.origin - p0;
        T a2 = a.mag2();
        T b2 = b.mag2();
        T d  = a.dot(b);
        T s  = d / a2;  // fractional projection of ray origin onto axis
        T x2 = d * s;
        T m  = ray.direction.dot(a);
        
        // construct quadratic coeffs
        T k0 = ray.direction.mag2() - m * m / a2;
        T k1 = 2 * (b.dot(ray.direction) - d * m / a2);
        T k2 = b2 - x2 - radius * radius;
        
        // setup for cylinder intersection
        T roots[2];
        Rect<T,1> interval;
        
        // construct an interval holding the ray's overlap with the infinite cylinder
        if (k0 != 0 and quadratic_solve(roots, k0, k1, k2)) {
            interval = Rect<T,1>::spanning_corners(roots[0], roots[1]);
        } else if (k2 < 0) {
            // if there is no intersection, then it might be the case that 
            // the ray points exactly down the interior of the cylinder, which
            // can only happen if the ray origin is inside too (k2 < 0). also
            // due to precision this manifests sometimes when there is one root:
            // when v and a are aligned then k0 is 0, and the equation can take the form
            // `0x^2 + (small)x + c = 0`; this should be treated the same as 'no solution'.
            // (Tangent rays with a single root will have k0 ≠ 0 and should be intersected
            // normally).
            interval = Rect<T,1>::full;
        } else {
            // no intersection
            return Rect<T,1>::empty;
        }
        
        // trace the endcaps
        // xxx: this is the wrong interval logic
        interval &= Sphere<T,N>(p0, radius).intersect(ray);
        interval &= Sphere<T,N>(p1, radius).intersect(ray);
        
        return interval;
    }
    */
    
};

} /* namespace geom */


template <typename T, index_t N>
struct std::hash<geom::Capsule<T,N>> {
    size_t operator()(const geom::Capsule<T,N> &c) const {
        constexpr size_t nonce = (size_t) 0x800065016103dd5fULL;
        return geom::hash_combine(
            geom::hash(c.p0),
            geom::hash(c.p1),
            geom::hash(c.radius)
        ) ^ nonce;
    }
};
