/*
 * LCRand.h
 *
 *  Created on: Nov 11, 2012
 *      Author: tbabb
 */

//TODO: test this

#ifndef LCRAND_H_
#define LCRAND_H_

#include <geomc/random/Random.h>

namespace geom {

/**
 * @ingroup random
 * @brief Implements a linear congruential pseudorandom number generator with 
 * period 2<sup>64</sup>. 
 * 
 * The parameters are the same as those chosen for Donald Knuth's MMIX language:
 * 
 *     multiplier: 6364136223846793005
 *     increment:  1442695040888963407
 * 
 * The underlying state is 64 bits; rand32() returns the 32 most significant digits.
 * 
 * The results of this generator are very fast and generally higher-quality
 * and longer-period than those returned by std::rand().
 * 
 * This generator is **NOT** cryptographically secure; do **NOT** use for 
 * gambling or cryptography. An attacker may predict all future numbers by 
 * observing one or two results.
 */
class LCRand: public geom::Random {
public:
    LCRand();
    LCRand(uint64_t seed);
    virtual ~LCRand();

    //returns 32 random bits.
    virtual uint32_t rand32();
    virtual void rseed(uint64_t seed);

    /// @return 64 random bits.
    inline uint64_t rand64() {
        uint64_t A = 6364136223846793005LL;
        uint64_t C = 1442695040888963407LL;
        state = A * state + C;
        return state;
    }
    
private:
    uint64_t state;
    static int uniq; //used to prevent successive unseeded constructions from
                     //mapping to the same state due to clock collision.
};

} /* namespace geom */
#endif /* LCRAND_H_ */
