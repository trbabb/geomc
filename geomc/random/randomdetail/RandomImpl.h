/* 
 * File:   RandomImpl.h
 * Author: tbabb
 *
 * Created on May 31, 2014, 6:07 PM
 */

#ifndef RANDOMIMPL_H
#define	RANDOMIMPL_H

namespace geom {
namespace detail {    
    
    template <typename T>
    class RandomImpl {
        // no implementation in the general case.
    };
    
    // partial specializations for rand<T>() must go here as partially
    // specialized RandomImpl implementations, since functions themselves 
    // cannot be partially specialized. thanks, c++.
    
}
}


#endif	/* RANDOMIMPL_H */

