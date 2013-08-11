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
 *  @brief Shape-related classes and functions.
 */

namespace geom {

    enum ArrayOrder {
        ARRAYORDER_FIRST_DIM_CONSECUTIVE, //colmajor; "fortran" order; if using (row,col) dimension order
        ARRAYORDER_LAST_DIM_CONSECUTIVE   //rowmajor; "C" order; if using (row, col) dimension order
    };
    
    template <typename T, index_t N> class Plane;
    template <typename T, index_t N> class Sphere;
    template <typename T, index_t N> class Bounded;
    template <typename T, index_t N> class Rect;
    template <typename T, index_t N> class Hit;
    
    template <typename T, index_t N, ArrayOrder Order=ARRAYORDER_FIRST_DIM_CONSECUTIVE> class GridIterator;
    
    typedef Sphere<double,2> Circle2d;
    typedef Sphere<float,2>  Circle2f;
    
    typedef Sphere<float,3>  Sphere3f;
    typedef Sphere<float,4>  Sphere4f;
    
    typedef Sphere<double,3> Sphere3d;
    typedef Sphere<double,4> Sphere4d;

    typedef Rect<double,1> Ranged;
    typedef Rect<double,1> Rect1d;
    typedef Rect<double,2> Rect2d;
    typedef Rect<double,3> Rect3d;
    typedef Rect<double,4> Rect4d;

    typedef Rect<index_t,1> Rangei;
    typedef Rect<index_t,1> Rect1i;
    typedef Rect<index_t,2> Rect2i;
    typedef Rect<index_t,3> Rect3i;
    typedef Rect<index_t,4> Rect4i;

    typedef Rect<float,1> Rangef;
    typedef Rect<float,1> Rect1f;
    typedef Rect<float,2> Rect2f;
    typedef Rect<float,3> Rect3f;
    typedef Rect<float,4> Rect4f;

};

#endif /* SHAPETYPES_H_ */
