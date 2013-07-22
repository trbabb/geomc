/*
 * Intersection.h
 *
 *  Created on: Mar 7, 2013
 *      Author: tbabb
 */

#ifndef INTERSECTION_H_
#define INTERSECTION_H_

namespace geom {

template <typename T, index_t N>
inline bool intersects(const Plane<T,N> &p, const Sphere<T,N> &s) const {
    return p.distance(s.center) >= s.radius;
}

template <typename T, index_t N>
inline bool intersects(const Sphere<T,N> &s, const Plane<T,N> &p) const {
    return intersects(p,s);
}

template <typename T, index_t N>
inline bool intersects(const Plane<T,N> &p, const Rect<T,N> &b) const {
    T dist = 0;
    for (index_t i = 0; i < (1 << N); i++){
        //TODO: this.
    }
    return false;
}

};

#endif /* INTERSECTION_H_ */
