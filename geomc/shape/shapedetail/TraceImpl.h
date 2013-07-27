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

template<typename T, index_t N>
bool _ImplTracePlane(T* s, const Plane<T,N> &p, const Ray<T,N> &r, HitSide *side){
    T dot = p.normal.dot(r.direction);
        
    if (dot < 0 and (*side & HIT_FRONT)){
        *side = HIT_FRONT;
    } else if (dot >= 0 and (*side & HIT_BACK)){
        *side = HIT_BACK;
    } else {
        // backface culled.
        return false;
    }
    
    T d = p.distance(r.origin); 
    *s = -d / dot; // similar triangles to get to the plane
    
    if (*s <= 0 || dot == 0){
        // hit point behind ray origin, or ray perfectly parallel to plane (prevent NaN)
        return false;
    }
    
    return true;
}

}; // namespace detail
}; // namespace geom


#endif /* SHAPEIMPL_H_ */
