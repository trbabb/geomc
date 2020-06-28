/*
 * Sphere.h
 *
 *  Created on: May 10, 2009
 *      Author: Tim Babb
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include <geomc/linalg/Vec.h>
#include <geomc/shape/shapedetail/Hit.h>
#include <geomc/shape/Bounded.h>
#include <geomc/function/Utils.h>

#include "shapedetail/TraceImpl.h"

namespace geom {

    /** 
     * @ingroup shape
     * @brief A N-dimensional circle, sphere, or hypersphere
     */
    template <typename T, index_t N>
    class Sphere : virtual public Convex<T,N> {
    public:
        /// Center of the sphere.
        Vec<T,N> center;
        /// Radius of the sphere.
        T r;

        /**
         * Construct a sphere at the origin with radius 1.
         */
        Sphere():center(Vec<T,N>::zeros),r(1) {}
        
        /**
         * Construct a sphere with center at the point `c`, having radius `r`.
         * @param c Center of sphere.
         * @param r Radius of spehre.
         */
        Sphere(Vec<T,N> c, T r):center(c),r(r) {}
        
        /**
         * @return An axis-aligned bounding box completely containing this
         * sphere.
         */
        Rect<T,N> bounds() {
            Vec<T,N> rvec = Vec<T,N>(r);
            return Rect<T,N>(center-rvec, center+rvec);
        }

        /**
         * Sphere-point intersection test.
         * @param p A point.
         * @return `true` if `p` is inside or on the surface of the sphere, `false`
         * otherwise.
         */
        bool contains(Vec<T,N> p) const {
            return center.dist2(p) <= r*r;
        }
        
        /**
         * Sphere-sphere intersection test.
         * @param s Another sphere.
         * @return `true` if `s` overlaps with this sphere's volume, false otherwise.
         */
        bool intersects(Sphere s) const {
            return s.center.dist2(center) <= r*r;
        }
        
        Vec<T,N> convex_support(Vec<T,N> d) const {
            return center + d.unit() * r;
        }
        
        /**
         * Sphere-ray intersection test.
         * @param ray The ray to intersect with this sphere.
         * @param side Whether to hit-test the front (outside) or back-facing (inside)
         * surfaces of this sphere.
         * @return A ray hit describing whether and where the ray intersects this sphere,
         * as well as the normal, side hit, and ray parameter.
         */
        Hit<T,N> trace(const Ray<T,N> &ray, HitSide side) const {
            Hit<T,N> h = Hit<T,N>(ray, side); // defaults to miss 
            T r2 = r*r;
            // if inside, we are guaranteed to hit the back. return miss if not tracing backface.
            if ((ray.origin.dist2(center) < r2) and !(side & HIT_BACK)) return h;
            const Vec<T,N> &dir = ray.direction; //for shorthand
            Vec<T,N> x0 = center - ray.origin;
            // solve for s such that ||s*ray - ctr|| == radius 
            T a = dir.dot(dir);
            T b = -2*dir.dot(x0);
            T c = x0.dot(x0) - r2;
            T roots[2];
            if (quadratic_solve(roots,a,b,c)) {
                T s;
                if (detail::chooseRayHit(&s, roots, &side)) {
                    // successful hit
                    h.p = ray.atMultiple(s);
                    h.n = (h.p - center).unit();
                    h.s = s;
                    h.side = side;
                    h.hit  = true;
                    return h;
                }
            } else {
                // no intersection; return miss.
                return h;
            }
        }
        
    }; /* Sphere<T,N> */
    
    
    /**
     * Intersect a ray and sphere, returning the ray multiples where the two shapes
     * intersect in `s0` and `s1`, or return `false` if there is no intersection.
     * @param ray Ray to intersect. `s0` or `s1` may be negative, and are ordered 
     * lowest to highest.
     * 
     * @param [out] s0 Lowest ray multiple generating a hit.
     * @param [out] s1 Highest ray multiple generating a hit.
     * @param center Center of sphere to trace.
     * @param r Radius of sphere to trace.
     * @return `true` if the ray intersects the sphere, `false` otherwise.
     */
    template <typename T, index_t N>
    bool trace_sphere(const Ray<T,N> &ray, T *s0, T *s1, const Vec<T,N> &center, T r) {
        T r2 = r*r;
        const Vec<T,N> &dir = ray.direction; //for shorthand
        Vec<T,N> x0 = center - ray.origin;
        // solve for s such that ||s*ray - ctr|| == radius 
        T a = dir.dot(dir);
        T b = -2 * dir.dot(x0);
        T c = x0.dot(x0) - r2;
        T roots[2];
        if (quadratic_solve(roots,a,b,c)) {
            // order the roots along the ray, from -inf; s0 first.
            if (roots[0] < roots[1]) {
                *s0 = roots[0];
                *s1 = roots[1];
            } else {
                *s0 = roots[1];
                *s1 = roots[0];
            }
            return true;
        }
        return false;
    }
    
} /* namespace geom */
#endif /* SPHERE_H_ */
