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

namespace geom {

template <typename T, index_t N>
class Plane {
public:
    
    // plane equation: ax + by + cz + d = 0
    // normal = (a,b,c)
    // d = distance along normal from origin to plane
    //     negative means move forward along normal,
    //     since we will have to go "more positive"
    //     to reach zero.
    
    Vec<T,N> normal; // unit normal
    T d; // distance along unit normal from origin to plane
         

    /*****************************
     * Structors                 *
     *****************************/
    
    Plane():d(0) { normal[0] = 1; }
    
    Plane(Vec<T,N> n):normal(n.unit()),d(0) {}
    
    Plane(Vec<T,N> n, Vec<T,N> origin):normal(n.unit()),d(-origin.dot(normal)) {}
    
    static inline Plane from_basis(Vec<T,N> a, Vec<T,N> b, 
            typename boost::enable_if_c<N == 3, int>::type = 0) {
        return Plane(a.cross(b));
    }
    
    static inline Plane from_basis(Vec<T,N> a, Vec<T,N> b, Vec<T,N> origin,
            typename boost::enable_if_c<N == 3, int>::type = 0) {
        return Plane(a.cross(b), origin);
    }
    
    static inline Plane from_triangle(Vec<T,N> p0, Vec<T,N> p1, Vec<T,N> p2,
            typename boost::enable_if_c<N == 3, int>::type = 0) {
        return Plane((p1-p0).cross(p2-p0), p0);
    }
    
    static inline Plane from_triangle(Vec<T,N> p[3], 
            typename boost::enable_if_c<N == 3, int>::type = 0) {
        return from_triangle(p[0], p[1], p[2]);
    }
    
    /*****************************
     * Methods                   *
     *****************************/
    
    // returns an s such that p - N*s is on the plane,
    // or what multiple of the normal must I travel off the plane?
    // in other words, positive "above" the plane (in the
    // direction of the normal), and negative "below".
    inline T distance(const Vec<T,N> &p) const {
        return normal.dot(p) + d;
    }

    inline bool contains(const Vec<T,N> &p) const {
        return distance(p) < 0;
    }
    
    inline Vec<T,N> project(const Vec<T,N> &p) const {
        return -(p.dot(normal) + d) * normal + p; 
    }
    
    // point on plane closest to origin
    inline Vec<T,N> origin() const {
        return -d * normal;
    }
    
    void apply(const AffineTransform<T,N> &xf) {
        Vec<T,N> p0 = normal*d; // a point on the plane
        Vec<T,N> p1 = xf.apply(p0); // the new position of that point
        normal = xf.applyNormal(normal).unit(); // construct a plane with transformed normal
        d = p1.dot(normal); // and transformed position
    }
    
    void applyInverse(const AffineTransform<T,N> &xf) {
        Vec<T,N> p0 = normal*d;
        Vec<T,N> p1 = xf.applyInverse(p0);
        normal = xf.applyInverseNormal(normal).unit();
        d = p1.dot(normal);
    }
    
    inline bool intersects(const Vec<T,N> *hull, index_t npts) const {
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
    
    inline bool intersects(const Plane<T,N> &p) const {
        return (p.normal.cross(normal) != 0);
    }
     
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
    
    friend Plane<T,N> operator*(const AffineTransform<T,N> &at, Plane<T,N> p) {
        p.apply(at);
        return p;
    }
    
    friend Plane<T,N> operator/(Plane<T,N> p, const AffineTransform<T,N> &at) {
        p.applyInverse(at);
        return p;
    }
    
    Plane<T,N>& operator*=(const AffineTransform<T,N> &tx) {
        this->apply(tx);
        return *this;
    }
    
    Plane<T,N>& operator/=(const AffineTransform<T,N> &tx) {
        this->applyInverse(tx);
        return *this;
    }
    
}; /* Plane<T,N> */

};

#endif /* PLANE_H_ */
