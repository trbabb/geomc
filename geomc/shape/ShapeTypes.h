/*
 * ShapeTypes.h
 *
 *  Created on: Nov 11, 2010
 *      Author: tbabb
 */

#ifndef SHAPETYPES_H_
#define SHAPETYPES_H_

#include <boost/utility/enable_if.hpp>
#include <geomc/linalg/LinalgTypes.h>

/** @defgroup shape Shape
 *  @brief Shape-related functions and classes.
 */

namespace geom {
    
    /** @addtogroup shape 
     *  @{
     */


    /** 
     * @brief Array traversal order, specified in terms of which axes to increment first.
     */
    enum ArrayOrder {
        /**
         * @brief Array traversal by incrementing the first coordinate in the 
         * innermost loop. 
         * 
         * For matrices, whose coordinates are ordered `(row, col)`, this 
         * represens column-major (i.e. "fortran") order. If coordinates
         * are `(x, y, ...)`, then this is row-major order.
         */
        ARRAYORDER_FIRST_DIM_CONSECUTIVE,
        /**
         * @brief Array traversal by incrementing the last coordinate in the
         * innermost loop.
         * 
         * For matrices, whose coordinates are ordered `(row, col)`, this
         * represents row-major (i.e. "C") order. If coordinates are `(x, y, ...)`,
         * then this is column-major order.
         */
        ARRAYORDER_LAST_DIM_CONSECUTIVE
    };
    
    // base shape concepts
    template <typename T, index_t N, typename Derived> class Bounded;
    template <typename T, index_t N, typename Derived> class Convex;
    template <typename T, index_t N, typename Derived> class Projectable;
    template <typename T, index_t N, typename Derived> class SdfEvaluable;
    template <typename T, index_t N, typename Derived> class RayIntersectable;
    
    // virtualization helper
    template <typename T, index_t N> class AnyConvex;
    template <typename Shape>        class AnyConvexImpl;
    
    // shapes
    template <typename T, index_t N> class Plane;
    template <typename T, index_t N> class Sphere;
    template <typename T, index_t N> class Rect;
    template <typename T, index_t N> class Cylinder;
    template <typename T, index_t N> class Simplex;
    template <typename Shape>        class Extrusion;
    template <typename Shape>        class Oriented;
    template <typename Shape>        class Frustum;
    template <typename Shape>        class Dilated;
    
    template <typename T, index_t N, typename Object, typename NodeData=void*>
        class KDTree;
    template <typename T, index_t N, ArrayOrder Order=ARRAYORDER_FIRST_DIM_CONSECUTIVE>
        class GridIterator;

/// @} // addtogroup shape
    
} // namespace geom

#endif /* SHAPETYPES_H_ */
