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

template <typename T, index_t N> class Bounded {
public:
    Bounded() {
        //do nothing
    }
    
    virtual ~Bounded(){
        //do nothing
    }

    virtual geom::Rect<T,N> bounds() = 0;
};

}
#endif /* BOUNDED_H_ */
