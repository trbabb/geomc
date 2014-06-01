/* 
 * File:   RandomDual.h
 * Author: tbabb
 * 
 * This class provides random number generation for duals when both
 * Dual.h and Random.h are in play.
 *
 * Created on December 23, 2013, 1:00 PM
 */

#ifndef RANDOMDUAL_H
#define	RANDOMDUAL_H

namespace geom {
namespace detail {
    
    template <typename T>
    class RandomImpl< Dual<T> > {
        
        public:
        
        static Dual<T> rand(Random *rng) {
            return Dual<T>(rng->rand<T>(), 0);
        }
        
        static Dual<T> rand(Random *rng, const Dual<T> &max) {
            return Dual<T>(rng->rand<T>(max.x), 0);
        }
        
        static Dual<T> rand(Random *rng, const Dual<T> &lo, const Dual<T> &hi) {
            return Dual<T>(rng->rand<T>(lo.x, hi.x), 0);
        }
        
    };
    
}
}


#endif	/* RANDOMDUAL_H */

