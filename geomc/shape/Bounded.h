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
template <typename T, index_t _N> class Bounded {
public:
    
    /// The coordinate type of this shape.
    typedef T elem_t;
    /// The dimension of this shape.
    static constexpr size_t N = _N;
    
    Bounded() {}
    virtual ~Bounded() {}

    /**
     * @return An axis-aligned box completely enclosing this shape.
     */
    virtual Rect<T,N> bounds() = 0;
};


/// Virtual class describing convex shapes in N-dimensional space.
template <typename T, index_t N> class Convex : virtual public Bounded<T,N> {
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
    virtual point_t convex_support(point_t d) const = 0;
    
    /**
     * @return True if and only if this convex shape overlaps `other`; false otherwise.
     * 
     * Note that some subclasses have overloaded versions of this function for specific
     * shapes which may offer better performance. When handling pointers to Convex objects,
     * it may be adviseable to cast to the derived class if it is known.
     * 
     * @param other Any other convex shape.
     */
    bool intersects(const Convex<T,N> &other) const {
        return gjk_intersect(*this, other);
    }
    
    /**
     * @brief Returns an axis-aligned box completely enclosing this shape.
     *
     * The default implementation calls `convex_support()` along each of the 
     * principal axes to find the extents.
     */
    Rect<T,N> bounds() {
        Rect<T,N> b;
        T* mins = PointType<T,N>::iterator(b.lo);
        T* maxs = PointType<T,N>::iterator(b.hi);
        for (index_t i = 0; i < N; ++i) {
            point_t axis, x;
            PointType<T,N>::iterator(axis)[i] =  1;
            x = convex_support(axis);
            maxs[i] = PointType<T,N>::iterator(x)[i];
            PointType<T,N>::iterator(axis)[i] = -1;
            x = convex_support(axis);
            mins[i] = PointType<T,N>::iterator(x)[i];
        }
        return b;
    }
    
};

/// @} // addtogroup shape

}
#endif /* BOUNDED_H_ */
