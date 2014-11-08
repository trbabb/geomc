/*
 * Bounded.h
 *
 * A <Bounded> is a thing which can bound its dimensions with an axis-aligned bounding box
 * via RectBound.
 *
 *  Created on: Oct 7, 2010
 *      Author: tbabb
 */

#ifndef BOUNDED_H_
#define BOUNDED_H_

#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Rect.h>

namespace geom {
    
/** @addtogroup shape
 *  @{
 */


/// Virtual class describing shapes with finite extents in N dimensions.
template <typename T, index_t N> class Bounded {
public:
    Bounded() {}
    virtual ~Bounded() {}

    /**
     * @return An axis-aligned box completely enclosing this shape.
     */
    virtual geom::Rect<T,N> bounds() = 0;
};


/// Virtual class describing convex shapes in N-dimensional space.
template <typename T, index_t N> class Convex {
public:
    typedef typename PointType<T,N>::point_t point_t;
    
    Convex() {}
    virtual ~Convex() {};
    
    /**
     * Returns the point on the surface of this convex shape that is furthest 
     * along direction `d` (i.e., has the highest dot product with `d`). 
     * 
     * All shapes which implement this function automatically support geometrical
     * intersection tests with any other Convex object.
     * 
     * @param d Direction along which to find a support plane.
     * @return A point on the surface of this convex shape.
     */
    virtual point_t convexSupport(point_t d) const = 0;
    
    /**
     * @return True if and only if this convex shape overlaps `other`; false otherwise.
     * 
     * Note that some subclasses have overloaded versions of this function for specific
     * shapes which may offer better performance. When handling points to Convex objects,
     * it may be adviseable to cast to the derived classes.
     * 
     * @param other Any other convex shape.
     */
    bool intersects(const Convex<T,N> &other) const {
        return gjk_intersect(*this, other);
    }
    
};

/// @} // addtogroup shape

}
#endif /* BOUNDED_H_ */
