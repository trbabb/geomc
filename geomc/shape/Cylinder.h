/* 
 * File:   Cylinder.h
 * Author: tbabb
 *
 * Created on May 19, 2014, 10:58 PM
 */

#ifndef CYLINDER_H
#define	CYLINDER_H

#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Plane.h>
#include <geomc/shape/Bounded.h>
#include <geomc/function/Utils.h>

#include "shapedetail/TraceImpl.h"

namespace geom {
    
    /** 
     * @ingroup shape
     * @brief An N-dimensional cylinder, given by its radius and endpoints.
     * 
     * Represents an extrusion of an N-2 sphere. In other words, an extrusion of 
     * a disk in 3D; an extrusion of a sphere in 4D; and an extrusion of a line 
     * segment in 2D.
     */
    template <typename T, index_t N>
    class Cylinder : virtual public Convex<T,N> {
        public:
            /// Axis endpoint.
            Vec<T,N> p0;
            /// Axis endpoint.
            Vec<T,N> p1;
            /// Cylinder radius.
            T radius;
            
            /// Construct a cylinder of radius and length 1, with axis along X+.
            Cylinder() : radius(1) { p1[0] = 1; }
            
            /// Construct a cylinder of radius `r` and length 1, with axis along X+.
            Cylinder(T r) : radius(r) { p1[0] = 1; }
            
            /**
             * @brief Construct a cylinder with arbitrary radius and endpoints.
             * @param p0 An endpoint of the cylinder axis.
             * @param p1 An endpoint of the cylinder axis.
             * @param radius Radius of cylinder.
             */
            Cylinder(const Vec<T,N> &p0, const Vec<T,N> &p1, T radius) :
                    p0(p0),
                    p1(p1),
                    radius(radius) {}
            
            /**
             * Cylinder-point intersection test.
             * @param p A point.
             * @return `true` if `p` is inside or on the surface of this cylinder;
             * `false` otherwise.
             */
            bool contains(const Vec<T,N> &p) const {
                Vec<T,N> v = p1 - p0;
                Vec<T,N> b = p - p0;
                T b2 = b.mag2();
                T v2 = v.mag2();
                T d  = v.dot(b);
                T x2 = d * d / v2;
                bool in_r = b2 - x2 <= radius * radius;
                bool in_caps = x2 < v2 && d >= 0;
                return in_r && in_caps;
            }
            
            /**
             * @return An axis-aligned bounding box completely containing this
             * cylinder.
             */
            Rect<T,N> bounds() {
                // construct two bounding boxes, one for each disk cap.
                // we do this by finding the extent of the cap along each axis.
                Vec<T,N> v  = p1 - p0;
                // square of unit normal's components
                // aka square of projection of normal onto axis of interest
                Vec<T,N> n2 = (v * v) / v.mag2(); 
                // solve for the remaining side of the triangle
                // which happens to be the projection of a perpendicular vector
                Vec<T,N> c0 = radius * std::sqrt((Vec<T,N>::ones - n2));
                // construct the disk bounds at the end cap positions
                Rect<T,N> b0(p0 - c0, p0 + c0);
                Rect<T,N> b1(p1 - c0, p1 + c0);
                // box union
                return b0 | b1;
            }
            
            Vec<T,N> convex_support(Vec<T,N> d) const {
                Vec<T,N> v = p1 - p0;
                Vec<T,N> perp = (d - d.projectOn(v)).unit() * radius;
                return (d.dot(v) > 0 ? p1 : p0) + perp;
            }
            
            /**
             * Ray-cylinder intersection.
             * @param ray A ray.
             * @param side Whether to hit-test the front (outside) or back 
             * (inside) of the cylinder.
             * @return A ray Hit describing whether and where the ray intersects 
             * this cylinder, as well as the normal, side hit, and ray parameter.
             */
            Hit<T,N> trace(const Ray<T,N> &ray, HitSide side) const {
                Hit<T,N> h;
                /* The general principle is as follows:
                 * We examine the right triangle formed by the candidate hit point P,
                 * `p0` and the cylinder axis. We may solve one of the triangle
                 * legs to obtain the distance from P to the axis. With P constructed
                 * as `o + sv`, we may solve for s such that the square of said
                 * distance equals the square of the cylinder radius. We then use CSG
                 * with two planes to remove the cylinder beyond the endcaps.
                 */
                
                Vec<T,N> axis = p1 - p0;
                Vec<T,N> k = ray.origin - p0;
                T l  = k.dot(axis);
                T m  = ray.direction.dot(axis);
                T v2 = ray.direction.mag2();
                T a2 = axis.mag2();
                T k2 = k.mag2();
                
                T a  = v2 - m * m / a2;
                T b  = 2 * (k.dot(ray.direction) - l * m / a2);
                T c  = k2 - l * l / a2 - radius * radius;
                
                T roots[2];
                if (quadratic_solve(roots, a, b, c)) {
                    // trace endcaps
                    T pl_0 = trace_plane(p0, -axis, ray);
                    T pl_1 = trace_plane(p1,  axis, ray);
                    T n_sign = 1;
                    
                    // order roots low -> hi
                    if (pl_0 > pl_1) {
                        std::swap(pl_0, pl_1); 
                        n_sign = -1;
                    }
                    if (roots[0] > roots[1]) {
                        std::swap(roots[0], roots[1]);
                    }
                    
                    // perform csg
                    bool hit_cylinder;
                    bool hit_near;
                    h.hit = csg_intersect(
                                roots[0], roots[1],
                                pl_0, pl_1,
                                &h.s,
                                &hit_cylinder,
                                &hit_near,
                                &h.side);
                    
                    // update the hit with P, N
                    // (csg_intersect took care of the others)
                    if (h.hit) {
                        h.p = ray.atMultiple(h.s);
                        if (hit_cylinder) {
                            h.n = Ray<T,N>(p0, axis).directionFromAxisTo(h.p).unit();
                        } else {
                            h.n = (((hit_near ? -1 : 1) * n_sign) * axis).unit();
                        }
                    }
                }
                
                return h;
            }
            
    };
    
    
} // namespace geom

#endif	/* CYLINDER_H */

