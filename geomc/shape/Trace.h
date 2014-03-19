/*
 * Trace.h
 *
 *  Created on: Mar 6, 2013
 *      Author: tbabb
 */

#ifndef TRACE_H_
#define TRACE_H_

#include <geomc/shape/shapedetail/Hit.h>
#include <geomc/shape/shapedetail/TraceImpl.h>
#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Plane.h>

namespace geom {


/********************************
 * Static tracing functions     *
 ********************************/

/*
template <typename T>
Hit<T,3> trace_bilinearPatch(Vec<T,3> p[4], const Ray<T,3> &r){
    Vec<T,3> a,b,c,d;
    a = p[3] - p[2] - p[1] + p[0];
    b = p[2] - p[0];
    c = p[1] - p[0];
    d = p[0];
    
}*/

//TODO: planar figures, to plug in with trace_planar_quad

template <typename T, index_t N>
inline Hit<T,3> trace_tri(Vec<T,3> p[3], const Ray<T,3> &r, HitSide sides=HIT_FRONT){
    return trace_tri(p[0], p[1], p[2], r, sides);
}

template <typename T>
Hit<T,3> trace_tri(const Vec<T,3> &p0, 
                   const Vec<T,3> &p1, 
                   const Vec<T,3> &p2, 
                   const Ray<T,3> &r, 
                   HitSide sides=HIT_FRONT){
    Hit<T,3> hit(r,sides);
    
    Vec<T,3> dPdu = p1 - p0;
    Vec<T,3> dPdv = p2 - p0;
    Plane<T,3> plane = Plane<T,3>::from_basis(dPdu, dPdv, p0);
    T s;
    
    if (!detail::_ImplTracePlane(&s, plane, r, &sides)) {
        return hit;
    }
    
    Vec<T,3> p = r.atMultiple(s);
    
    // dPdu, dPdv, and N form a 3d basis.
    // thus we can solve:
    //    u*dPdu + v*dPdv + t*N = P
    // t should be very close to zero; we can ignore it.
    
    // 3d matrix inversion is faster than you think. 
    // also, gcc tends to unroll those copy ops.
    SimpleMatrix<T,3,3> mtx;
    std::copy(dPdu.begin(), dPdu.end(), mtx.col(0));
    std::copy(dPdv.begin(), dPdv.end(), mtx.col(1));
    std::copy(plane.normal.begin(), plane.normal.end(), mtx.col(2));
    inv(&mtx, mtx);
    Vec<T,3> uvn = mtx * (p-p0);
    
    T u = uvn.x;
    T v = uvn.y;
    
    if (u < 0 || v < 0 || u + v > 1){
         //intersection outside of triangle
        return hit;
    }
    
    hit.p = p;
    hit.n = plane.normal;
    hit.s = s;
    hit.side = sides; //set by ImplTracePlane
    hit.hit = true;
    
    return hit;
}

template <typename T>
Hit<T,3> trace_planar_quad(const Vec<T,3> &dPdu, const Vec<T,3> &dPdv, const Vec<T,3> &origin, const Ray<T,3> &ray, HitSide sides){
    Hit<T,3> hit(ray,sides);
    Plane<T,3> pl = Plane<T,3>::from_basis(dPdu, dPdv, origin);
    T s;
    
    if (!detail::_ImplTracePlane(&s, pl, ray, &sides)){
        // miss
        return hit;
    }
    
    Vec<T,3> p = ray.atMultiple(s);
    
    // solve for u and v
    SimpleMatrix<T,3,3> mtx;
    std::copy(dPdu.begin(), dPdu.end(), mtx.col(0));
    std::copy(dPdv.begin(), dPdv.end(), mtx.col(1));
    std::copy(pl.normal.begin(), pl.normal.end(), mtx.col(2));
    inv(&mtx, mtx);
    Vec<T,3> uvn = mtx * (p-origin);
    
    T u = uvn.x;
    T v = uvn.y;
    
    if (u < 0 or u > 1 or v < 0 or v > 1){
        // missed quad area
        return hit;
    }
    
    hit.p = p;
    hit.n = pl.normal;
    hit.s = s;
    hit.side = sides; // set by ImplTracePlane
    hit.hit = true;
    
    return hit;
}


} // namespace geom

#endif /* TRACE_H_ */
