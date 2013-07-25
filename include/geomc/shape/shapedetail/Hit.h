/*
 * Hit.h
 *
 *  Created on: Mar 8, 2013
 *      Author: tbabb
 */

#ifndef HIT_H_
#define HIT_H_

#include <limits>
#include <boost/utility/enable_if.hpp>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Ray.h>

namespace geom {

// there are several things we may independently want to query with a ray:
//   hit? yes/no
//   ray s
//   surface P
//   surface N
//   surface U and V
//   surface derivatives (dPdu, dPdv, ...)
//   world-to-object space transform of hit object.
// perhaps these should be exposed via overloaded trace() functions
// with different kinds of output Hit parameters.

// bitwise/masking logic allowed
// i.e. HIT_BOTH = HIT_FRONT | HIT_BACK

enum HitSide {
    HIT_FRONT = 1,
    HIT_BACK  = 2,
    HIT_BOTH  = 3
};

// TODO: allow mechanism for supplying other primvars.
//       (I think only U,V(,W?) are needed; plus surface derivatives (a matrix; dPdu, dPdv, dPdw);
//       other user-supplied primvars will be interped based on these anyway).
//       Note: W could just be N... maybe. for volume tracing, mebbe not.

template <typename T, index_t N>
struct Hit {
    Vec<T,N> p;
    Vec<T,N> n;
    T s;
    Ray<T,N> r;
    bool hit;
    HitSide side;
    
    Hit():s(std::numeric_limits<T>::max()),hit(false),side(HIT_BOTH){}
    
    Hit(Ray<T,N> r, HitSide s):s(std::numeric_limits<T>::max()),r(r),hit(false),side(s){}
};

}; // namespace geom

#endif /* HIT_H_ */
