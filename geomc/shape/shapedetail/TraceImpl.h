/*
 * ShapeImpl.h
 *
 *  Created on: Mar 8, 2013
 *      Author: tbabb
 */

#ifndef SHAPEIMPL_H_
#define SHAPEIMPL_H_

#include <geomc/shape/shapedetail/Hit.h>
#include <geomc/shape/ShapeTypes.h>

namespace geom {
namespace detail {

template <typename T, index_t N>
bool _ImplTracePlane(T* s, const Plane<T,N> &p, const Ray<T,N> &r, HitSide *side) {
    T dot = p.normal.dot(r.direction);
        
    if (dot < 0 and (*side & HIT_FRONT)) {
        *side = HIT_FRONT;
    } else if (dot >= 0 and (*side & HIT_BACK)) {
        *side = HIT_BACK;
    } else {
        // backface culled.
        return false;
    }
    
    T d = p.distance(r.origin); 
    *s = -d / dot; // similar triangles to get to the plane
    
    if (*s <= 0 || dot == 0) {
        // hit point behind ray origin, or ray perfectly parallel to plane (prevent NaN)
        return false;
    }
    
    return true;
}

// choose a ray hit for a fully convex object (i.e. one with only a single ray 
// entry and exit each), accounting for hit side preference.
template <typename T>
bool chooseRayHit(T *s, T roots[2], HitSide *side) {
    T s0, s1;
    
    // order the roots along the ray, from -inf; s0 first.
    if (roots[0] < roots[1]) {
        s0 = roots[0];
        s1 = roots[1];
    } else {
        s0 = roots[1];
        s1 = roots[0];
    }

    // lesser root will always be the frontside, greater will be the back.
    // if tracing both front and back, and both hit, return the nearer (front).
    if (s0 > 0 && (*side & HIT_FRONT)) {
        *s = s0;
        *side = HIT_FRONT;
    } else if (s1 > 0 && (*side & HIT_BACK)) {
        *s = s1;
        *side = HIT_BACK;
    } else {
        return false; // hit behind origin; return miss.
    }
    return true;
}

} // namespace detail
} // namespace geom

#endif /* SHAPEIMPL_H_ */
