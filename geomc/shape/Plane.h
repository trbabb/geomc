/*
 * Plane.h
 *
 *  Created on: Mar 7, 2013
 *      Author: tbabb
 * 
 * TODO: test affine transform ops
 * TODO: does AffineTransform * Vec<T,N+1>(normal, d) work?
 */

#ifndef PLANE_H_
#define PLANE_H_

#include <geomc/linalg/AffineTransform.h>
#include <geomc/shape/Trace.h>
#include <geomc/linalg/Orthogonal.h>

namespace geom {

/**
 * @ingroup shape
 * @brief A geometric plane or hyperplane
 * 
 * A Euclidean subspace of R<sup>N</sup> with dimension N-1.
 */
template <typename T, index_t N>
class Plane {
public:
    
    // plane equation: ax + by + cz + d = 0
    // normal = (a,b,c)
    // d = distance along normal from origin to plane
    //     negative means move forward along normal,
    //     since we will have to go "more positive"
    //     to reach zero.
    
    /// Normal direction, unit length.
    Vec<T,N> normal;
    /// Distance along the unit normal from the origin to the plane surface
    T d;
         

    /*****************************
     * Structors                 *
     *****************************/
    
    /**
     * Construct a plane with a normal along +x, crossing the origin.
     */
    Plane():d(0) { normal[0] = 1; }
    
    /**
     * Construct a plane with normal `n`.
     * @param n Normal direction. 
     */
    Plane(Vec<T,N> n):normal(n.unit()),d(0) {}
    
    /**
     * Construct a plane with normal `n` passing through the point `pt`.
     * @param n Normal direction.
     * @param pt A point on the plane.
     */
    Plane(Vec<T,N> n, Vec<T,N> pt):normal(n.unit()),d(-pt.dot(normal)) {}
    
    /**
     * Construct a plane spanning the given basis vectors. 
     * @param bases An array of `N-1` bases.
     * @return A plane whose normal is orthogonal to all the given bases.
     */
    static inline Plane from_basis(const Vec<T,N> bases[N-1]) {
        return Plane(orthogonal(bases));
    }
    
    /**
     * Construct a plane spanning the given basis vectors through the given point.
     * @param bases An array of `N-1` bases.
     * @param p Point through which to construct the plane.
     * @return A plane whose normal is orthogonal to all the given bases.
     */
    static inline Plane from_basis(const Vec<T,N> bases[N-1], Vec<T,N> p) {
        return Plane(orthogonal(bases), p);
    }
    
    /**
     * Construct a plane spanning the points on the given simplex (i.e. line 
     * segment, triangle, tetrahedron, etc; as appropriate for the dimension). 
     * @param points The `N` vertecies of a simplex.
     * @return A plane on which the given simplex lies.
     */
    static inline Plane from_simplex(const Vec<T,N> points[N]) {
        Vec<T,N> bases[N-1];
        for (index_t i = 1; i < N; i++) {
            // todo: can we choose the "pivot" more intelligently? 
            // i.e. avoid "sliver" simplexes.
            bases[i] = points[i] - points[0];
        }
        return Plane(bases, points[0]);
    }
    
    /*****************************
     * Methods                   *
     *****************************/
    
    /**
     * @return The distance from `p` along the normal direction to the surface
     * of the plane. Negative if `p` is on the backfacing side of the plane.
     */
    inline T distance(const Vec<T,N> &p) const {
        return normal.dot(p) + d;
    }

    /**
     * Half-space test.
     * @return `true` if `p` is on or below the surface of the plane (i.e. on
     * the side opposite the normal); `false` otherwise.
     */
    inline bool contains(const Vec<T,N> &p) const {
        return distance(p) <= 0;
    }
    
    /**
     * Half-space test.
     * @return `true` if `shape` is completely on or below the surface of the plane 
     * (i.e. on the side opposite the normal); false otherwise.
     */
    inline bool contains(const Convex<T,N> &shape) const {
        return shape.convexSupport(normal).dot(normal) <= -d;
    }
    
    /**
     * @return `p` projected onto the plane.
     */
    inline Vec<T,N> project(const Vec<T,N> &p) const {
        return -(p.dot(normal) + d) * normal + p; 
    }
    
    /**
     * @return The point on the plane closest to the origin.
     */
    inline Vec<T,N> origin() const {
        return -d * normal;
    }
    
    /**
     * Geometrically transform this plane by `xf`.
     * @param xf An affine transformation.
     */
    void apply(const AffineTransform<T,N> &xf) {
        Vec<T,N> p0 = normal*d; // a point on the plane
        Vec<T,N> p1 = xf.apply(p0); // the new position of that point
        normal = xf.applyNormal(normal).unit(); // construct a plane with transformed normal
        d = p1.dot(normal); // and transformed position
    }
    
    /**
     * Un-transform this plane by `xf`.
     * @param xf An affine transformation.
     */
    void applyInverse(const AffineTransform<T,N> &xf) {
        Vec<T,N> p0 = normal*d;
        Vec<T,N> p1 = xf.applyInverse(p0);
        normal = xf.applyInverseNormal(normal).unit();
        d = p1.dot(normal);
    }
    
    /**
     * Convex shape intersection test.
     * @return `true` if and only if this plane passes through `shape`.
     */
    inline bool intersects(const Convex<T,N> &shape) const {
        return shape.convexSupport( normal).dot(normal) + d  > 0 and
               shape.convexSupport(-normal).dot(normal) + d <= 0;
    }
    
    /**
     * Sphere interesection test.
     * @return `true` if and only if this plane passes through `s`.
     */
    inline bool intersects(const Sphere<T,N> &s) const {
        return distance(s.center) >= s.r;
    }
    
    /**
     * Convex hull intersection test.
     * @param hull A list of points whose convex hull is to be tested.
     * @param npts Number of points in the list.
     * @return `true` if this plane intersects the convex hull of `hull`, `false`
     * otherwise,
     */
    bool intersects(const Vec<T,N> *hull, index_t npts) const {
        bool negative;
        for (index_t i = 0; i < npts; i++){
            T dist = distance(hull[i]);
            if (i > 0){
                // if a point lies on the plane, or another point lies on the
                // other side of the plane, then this shape crosses the plane.
                if (dist == 0 || ((dist < 0) != negative)){
                    return false;
                }
            } else {
                negative = dist < 0;
            }
        }
        return false;
    }
     
    /**
     * Ray-plane intersection test.
     * @param ray A ray.
     * @param sides Whether to hit-test the front or back sides of this plane.
     * @return A ray Hit describing whether and where the ray intersects this Plane,
     * as well as the normal, side hit, and ray parameter.
     */
    Hit<T,N> trace(const Ray<T,N> &ray, HitSide sides) const {
        Hit<T,N> hit(ray, sides);
        T s;
        if (!detail::_ImplTracePlane(&s, *this, ray, &sides)) {
            // miss
            return hit;
        }
        hit.p = ray.atMultiple(s); 
        hit.n = normal;
        hit.s = s;
        hit.side = sides; // set by ImplTracePlane
        hit.hit = true;
        
        return hit;
    }
    
    /*****************************
     * Operators                 *
     *****************************/
    
    /**
     * Affine transformation.
     * @param tx An affine transformation.
     * @param p A plane.
     * @return A new plane representing `p` transformed by `at`.
     */
    friend Plane<T,N> operator*(const AffineTransform<T,N> &tx, Plane<T,N> p) {
        Plane<T,N> ret = p;
        ret.apply(tx);
        return ret;
    }
    
    /**
     * Inverse affine transformation.
     * @param p A plane.
     * @param tx An affine transformation.
     * @return A new plane representing `p` un-transformed by `at`.
     */
    friend Plane<T,N> operator/(Plane<T,N> p, const AffineTransform<T,N> &tx) {
        Plane<T,N> ret = p;
        ret.applyInverse(tx);
        return ret;
    }
    
    /**
     * Affine transformation.
     * 
     * @param tx An affine transformation to apply to this plane.
     */
    Plane<T,N>& operator*=(const AffineTransform<T,N> &tx) {
        this->apply(tx);
        return *this;
    }
    
    /**
     * Inverse affine transformation.
     * 
     * @param tx An affine transformation to un-apply to this plane.
     */
    Plane<T,N>& operator/=(const AffineTransform<T,N> &tx) {
        this->applyInverse(tx);
        return *this;
    }
    
}; /* Plane<T,N> */

} // namespace geom

#endif /* PLANE_H_ */
