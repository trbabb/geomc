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
    // no, the Hit should actually be called HitTest.
    // and the ray test function should be passed flags for what information
    // is queried. Then the Hit is written to accordingly. trace functions return
    // a bool.

// bitwise/masking logic allowed
// i.e. HIT_BOTH = HIT_FRONT | HIT_BACK

/** @ingroup shape 
 *  @{
 */
    
/** @brief Ray hit testing flags
 * 
 * May be bitwise-combined (i.e. `HIT_BOTH = HIT_FRONT | HIT_BACK`)
 */
enum HitSide {
    /// Ray-test front-facing surfaces
    HIT_FRONT = 1,
    /// Ray-test back-facing surfaces
    HIT_BACK  = 2,
    /// Ray-test front- and back-facing surfaces.
    HIT_BOTH  = 3
};

// TODO: allow mechanism for supplying other primvars.
//       (I think only U,V(,W?) are needed; plus surface derivatives (a matrix; dPdu, dPdv, dPdw);
//       other user-supplied primvars will be interped based on these anyway).
//       Note: W could just be N... maybe. for volume tracing, mebbe not.

/**
 * @ingroup shape 
 * @brief Ray hit class.
 */
template <typename T, index_t N>
struct Hit {
    /// Point of hit
    Vec<T,N> p;
    /// Normal of surface at hit point
    Vec<T,N> n;
    /// Ray parameter of hit, such that `r * s = p`
    T s;
    /// Ray which generated the hit
    Ray<T,N> r;
    /// Whether the ray hit the tested geometry
    bool hit;
    /// The side of the geometry which generated the hit.
    HitSide side;
    
    /// Construct a new ray miss.
    Hit():s(std::numeric_limits<T>::max()),hit(false),side(HIT_BOTH){}
    
    /// Construct a new ray miss with ray `r` and HitSides `s`.
    Hit(Ray<T,N> r, HitSide s):s(std::numeric_limits<T>::max()),r(r),hit(false),side(s){}
};

/// @} // ingroup linalg

} // namespace geom

#endif /* HIT_H_ */
