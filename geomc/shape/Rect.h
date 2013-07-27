/*
 * Rect.h
 *
 * An N-dimensional axis-aligned box.
 *
 * For N > 1, point_t is Vec<T,N>. For N=1, point_t is T;
 * thus a Rect<T,1> may be used to naturally represent a 1D range.
 *
 * Example: Rect<double,4>::point_t  is  Vec<double,4>;  getMin() returns Vec4d.
 *          Rect<double,1>::point_t  is  double;         getMin() returns a double.
 *
 * If you wanted to experiment with interesting types of bound (like T=Expr, double[], function, ...),
 * then you could use Rect<MyWeirdType, 1> (this is equivalent to the original RectBound implementation).
 *
 *  Created on: Oct 7, 2010
 *      Author: Tim Babb
 */

#ifndef Rect_H_
#define Rect_H_

#include <algorithm>
#include <cmath>
#include <functional>

#include <geomc/shape/Bounded.h>
#include <geomc/linalg/Ray.h>
#include <geomc/linalg/Vec.h>
#include <geomc/shape/shapedetail/Hit.h>

namespace geom {

/****************************
 * Rect Class               *
 ****************************/

template <typename T, index_t N> class Rect : virtual public Bounded<T,N> {
protected:
    typedef PointType<T,N> ptype;

public:
    typedef typename ptype::point_t point_t;

protected:
    point_t mins;
    point_t maxs;

public:

    Rect(point_t c1, point_t c2):
        mins(std::min(c1,c1)),
        maxs(std::max(c1,c2)) {
        //do nothing else
    }

    Rect():
        mins((T)0),
        maxs((T)1){
        //do nothing else
    }

    /****************************
     * Convenience functions    *
     ****************************/

    static Rect<T,N> fromCenter(typename Rect<T,N>::point_t center,
                                typename Rect<T,N>::point_t dims){
        return Rect<T,N>(center-dims/2, center+dims/2);
    }

    static Rect<T,N> fromCorner(typename Rect<T,N>::point_t corner,
                                typename Rect<T,N>::point_t dims){
        return Rect<T,N>(corner, corner+dims);
    }

    static bool contains(typename Rect<T,N>::point_t min,
                         typename Rect<T,N>::point_t max,
                         typename Rect<T,N>::point_t pt){
        min = std::min(min,max);
        max = std::max(min,max);
        for (index_t axis = 0; axis < N; axis++){
            T v = ptype::iterator(pt)[axis];
            if (v < ptype::iterator(min)[axis] or v >= ptype::iterator(max)[axis]){
                return false;
            }
        }
        return true;
    }

    /*****************************
     * Operators                   *
     *****************************/

    const Rect<T,N> operator+(const Rect<T,N> &b) const{
        return boxUnion(b);
    }

    const Rect<T,N> operator+(point_t dx) const{
        //translation
        return Rect<T,N>(mins + dx, maxs + dx);
    }

    Rect<T,N>& operator+=(const Rect<T,N> &b){
        //box union
        maxs = std::max(b.maxs, maxs);
        mins = std::min(b.mins, mins);
        return *this;
    }

    Rect<T,N>& operator+=(T dx){
        //translation
        maxs += dx;
        mins += dx;
        return *this;
    }

    bool operator==(Rect<T,N> b) const {
        return (maxs == b.maxs) && (mins == b.mins);
    }
    
    bool operator!=(Rect<T,N> b) const {
        return (maxs != b.maxs) || (mins != b.mins);
    }

    //uniform scale
    Rect<T,N> operator*(point_t a) const {
        return Rect<T,N>(mins*a, maxs*a);
    }

    //uniform scale
    Rect<T,N>& operator*=(point_t a) {
        mins *= a;
        maxs *= a;
        return *this;
    }

    Rect<T,N> operator/(point_t a) const {
        return Rect<T,N>(mins/a, maxs/a);
    }

    Rect<T,N>& operator/=(point_t a) {
        mins /= a;
        maxs /= a;
        return *this;
    }

    template <typename U, index_t M> operator Rect<U,M>(){
        return Rect<U,M>((typename Rect<U,M>::point_t) mins,
                         (typename Rect<U,M>::point_t) maxs);
    }


    /*****************************
     * Public Methods              *
     *****************************/

    bool contains(T pt) const {
        for (index_t axis = 0; axis < N; axis++){
            T v = ptype::iterator(pt)[axis];
            if (v < ptype::iterator(min)[axis] or v >= ptype::iterator(max)[axis]){
                return false;
            }
        }
        return true;
    }
 
    bool intersects(const Rect<T,N> &box) const {
        for (index_t axis = 0; axis < N; axis++){
            // disjoint on this axis?
            if (ptype::iterator(maxs)[axis]     <= ptype::iterator(box.mins)[axis] or 
                ptype::iterator(box.maxs)[axis] <= ptype::iterator(mins)[axis]){
                return false;
            }
        }
        return true;
    }

    inline point_t min() const {
        return mins;
    }

    inline point_t max() const {
        return maxs;
    }

    point_t getCenter() const {
        return (maxs+mins)/2;
    }

    inline point_t getDimensions() const {
        return maxs - mins;
    }

    void setDimensions(point_t dim) {
        dim = std::abs(dim);
        point_t diff = (dim - (maxs-mins)) / 2;
        mins = mins-diff;
        maxs = maxs+diff;
    }

    void setCorners(point_t corner1, point_t corner2) {
        maxs = std::max(corner1, corner2);
        mins = std::min(corner1, corner2);
    }
    
    void setCenter(point_t center) {
        point_t tx = center - getCenter();
        maxs += tx;
        mins += tx;
    }

    void translate(point_t tx) {
        maxs += tx;
        mins += tx;
    }
    
    T volume() const {
        T vol = 1;
        point_t dim = getDimensions();
        for (index_t axis = 0; axis < N; axis++){
            vol *= std::max(ptype::iterator(dim)[axis], 0);
        }
        return vol;
    }
    
    bool isEmpty() const {
        for (index_t axis = 0; axis < N; axis++){
            T hi = ptype::iterator(maxs)[axis];
            T lo = ptype::iterator(mins)[axis];
            if (hi <= lo){
                return true;
            }
        }
        return false;
    }
    
    Rect<T,N> rangeIntersection(const Rect<T,N> &b) const {
        return Rect<T,N>(std::max(mins, b.mins),
                         std::min(maxs, b.maxs));
    }

    Rect<T,N> rangeUnion(const Rect<T,N> &b) const {
        return Rect<T,N>(std::min(mins, b.mins),
                         std::max(maxs, b.maxs));
    }

    Rect<T,N> ptUnion(point_t pt) const {
        return Rect<T,N>(std::min(mins,pt),
                         std::max(maxs,pt));
    }

    //t is a coordinate on (0,1) representing a position within this box
    point_t boxCoords(point_t t) const {
        return (point_t(1)-t) * mins + t * maxs ;
    }
    
    Hit<T,N> trace(const Ray<T,N> &r, HitSide sides) const {
        T near_root = std::numeric_limits<T>::lowest(); // non-denormal extremes, for speed (as opposed to inf).
        T far_root  = std::numeric_limits<T>::max();    // besides, infinity may not be defined for T.
        index_t near_axis = 0;
        index_t far_axis  = 0;
        T near_coord  = 0; // for guaranteeing that the hit point is exactly on the rect surface,
        T far_coord   = 0; // regardless of precision and rounding error in the ray arithmetic.
        Hit<T,N> hit  = Hit<T,N>(r, sides); // defaults to miss
        
        // test ray intersection with an infinite slab along each axis
        // box intersections are the farthest near hit and nearest far hit
        for (index_t axis = 0; axis < N; axis++){
            if (r.direction[axis] == 0){
                // ray direction tangent to test planes, no intersection along this axis
                if (r.origin[axis] < mins[axis] || r.origin[axis] > maxs[axis]) {
                    // origin outside of test slab; miss
                    return hit;
                }
            } else {
                // coordinate of hit, along tested axis
                T c1 = maxs[axis];
                T c2 = mins[axis];
                // ray multiple of hit
                T s1 = (c1 - r.origin[axis]) / r.direction[axis];
                T s2 = (c2 - r.origin[axis]) / r.direction[axis];
                
                // order s1, s2 to (near, far) (aka front, back)
                if (s2 < s1){
                    std::swap(s1,s2);
                    std::swap(c1,c2);
                }
                if (s1 > near_root){
                    near_root  = s1;
                    near_axis  = axis; 
                    near_coord = c1;
                }
                if (s2 < far_root){
                    far_root  = s2;
                    far_axis  = axis;
                    far_coord = c2;
                }
            }
        }
        
        if (near_root > far_root){
            // miss
            return hit;
        } else if (near_root > 0 and (sides & HIT_FRONT)){
            // hit front
            hit.p = r * near_root;
            hit.p[near_axis] = near_coord;
            hit.s = near_root;
            hit.n[near_axis] = -copysign(1, r.direction[near_axis]);
            hit.side = HIT_FRONT;
            hit.hit = true;
        } else if (far_root > 0 and (sides & HIT_BACK)) {
            // hit back
            hit.p = r * far_root;
            hit.p[far_axis] = far_coord;
            hit.s = far_root;
            hit.n[far_axis] = copysign(1, r.direction[far_axis]);
            hit.side = HIT_BACK;
            hit.hit = true;
        }
        // if encountered cases above; hit.
        // else miss (box entirely behind ray origin)
        return hit;
    }

    /*****************************
     * Inherited Methods         *
     *****************************/

    Rect<T,N> bounds() {
        return *this;
    }

}; // end rect class

} //end namespace geom

//allows RectBound<T,N>s to work with std::tr1 hash containers
namespace std { namespace tr1 {

   //TODO: This probably isn't a great hashing function
   template <typename T, index_t N>
   struct hash< geom::Rect<T,N> > : public unary_function<geom::Rect<T,N>, size_t> {
       size_t operator()(const geom::Rect<T,N>& box) const {
           size_t h1 = std::tr1::hash<typename geom::Rect<T,N>::point_t>(box.getMin());
           size_t h2 = std::tr1::hash<typename geom::Rect<T,N>::point_t>(box.getMax());
           return h1 + h2 + (h1 * h2) + (h1 ^ h2);
       }
   };

}}

#endif /* Rect_H_ */
