/*
 * Sphere.h
 *
 *  Created on: May 10, 2009
 *      Author: Tim Babb
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "linalg/Vec.h"
#include "shapedetail/Hit.h"
#include "Bounded.h"
#include "Utils.h"

namespace geom {

    template <typename T, index_t N>
    class Sphere : virtual public Bounded<T,N> {
    public:
        Vec<T,N> center;
        T r;

        //structors
        Sphere():center(Vec<T,N>::zeros),r(1){}
        Sphere(Vec<T,N> c, T r):center(c),r(r){}
        
        Rect<T,N> bounds(){
            Vec<T,N> rvec = Vec<T,N>(r);
            return Rect<T,N>(center-rvec, center+rvec);
        }

        bool contains(Vec<T,N> p) const {
            return center.dist2(p) <= r*r;
        }
        
        bool intersects(Sphere s) const {
            return s.center.dist2(center) <= r*r;
        }
        
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
            if (quadratic_solve(roots,a,b,c)){
                T s, s0, s1;
                // order the roots along the ray, from -inf; s0 first.
                if (roots[0] < roots[1]){
                    s0 = roots[0];
                    s1 = roots[1];
                } else {
                    s0 = roots[1];
                    s1 = roots[0];
                }
                // lesser root will always be the frontside, greater will be the back.
                // if tracing both front and back and both hit, return the nearer (front).
                if (s0 > 0 && (side & HIT_FRONT)){
                    s = s0;
                    h.side = HIT_FRONT;
                } else if (s1 > 0 && (side & HIT_BACK)){
                    s = s1;
                    h.side = HIT_BACK;
                } else {
                    return h; // hit behind origin; return miss.
                }
                
                // successful hit
                h.p = ray.atMultiple(s);
                h.n = (h.p - center).unit();
                h.s = s;
                h.hit = true;
                return h;
            } else {
                // no intersection; return miss.
                return h;
            }
        }
        
    }; /* Sphere<T,N> */
    
}; /* namespace geom */
#endif /* SPHERE_H_ */
