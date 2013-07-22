/*
 * MersenneRand.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 * 
 * Mersenne twister pseudorandom number generator.
 * 
 * Advantages:
 *   - Very long period (does not repeat after a practically-attainable number of iterations)
 *   - High-quality random numbers (no patterns; matches uniform distribution)
 *   - Fast
 *   
 * Disadvantages:
 *   - NOT cryptographically secure!
 *     - Don't use for security/cryptography
 *     - Don't use for gambling
 */

//TODO: Test on random sphere; there seems to be "banding".
//TODO: Supply a way to fill the entire buffer with supplied entropy;
//      right now there are only 2^64 possible sequences. 

#include "Random.h"

#ifndef MERSENNERAND_H_
#define MERSENNERAND_H_

namespace geom {

    #define MT_LEN  624

    class MTRand: public Random {
    public:
        MTRand();
        MTRand(uint64_t seed);
        virtual ~MTRand();

        virtual void rseed(uint64_t seed);
        virtual uint32_t rand32();

    private:
        int mt_index;
        uint32_t mt_buffer[MT_LEN];
        
        static int uniq;
    };

}
#endif /* MERSENNERAND_H_ */
