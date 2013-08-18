/*
 * Random.h
 * 
 * (!) This class was designed to be portable between systems with different integer sizes,
 * although this feature has not been tested (!)
 * 
 *  Created on: Feb 22, 2009
 *  Heavily refactored in April 2011
 *      Author: Tim Babb
 */

//TODO: create a 3d random test suite.
//TODO: create an add_entropy() function, which hashes the bits and mixes them into the generator state.
//      Also be able to seed via the same methods as add_entropy.
//TODO: implement Yarrow.
//TODO: implement rand<char>(), and move 'uint32_t reserve_bits' to base class. implement getBits(char *dest, int bits)
//TODO: (low-pri) refactor rand<float> and rand<double> to use numerical limits for systems with half floats, etc. 

/**
 * @defgroup random
 * @brief Random number and object generation..
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <stdint.h>

namespace geom {
    
    /**
     * @ingroup random
     * @brief Abstraction for a source of random or pseudorandom bits.
     * 
     * Implementations provide random bits by overriding `rand32()`, which should 
     * return no fewer than 32 random bits.
     * 
     * Client-friendly random numbers are exposed via `Random::rand<T>()`. Specializations
     * are provided for a selection of primitive `T`. Specializations of `rand()` for 
     * other `T` may be written, which will automatically be inherited by all subclasses of `Random`.
     * 
     * The general contract of `rand()` which should be followed by all specializations:
     * 
     * `rand(void)`
     *   - If `T` is discrete, than this function should return values evenly over the entire 
     *     space of possible `T`.
     *   - If `T` is a continuous type, then this function should return values evenly
     *     over an analogue of the interval [0,1], inclusive.
     * 
     * `rand(T max)`
     *   - Should return random values evenly distributed over all possible states between
     *     0 (or the `T` analogue of 0) and the value of `max`
     * 
     * `rand(T lo, T hi)`
     *   - Should return random values evenly distributed over all possible states 
     *     bounded by `lo` and `hi`
     */
    class Random {
    public:

        /// Construct a new random number generator with an implementation-selected seed.
        Random(void);
        
        /// Construct a new random number generator with the bits of `seed` as a source
        /// of entropy.
        Random(uint64_t seed); //allow us to supply up to 64 bits of entropy
        virtual ~Random();

        /// @return No fewer than 32 (pseudo-) random bits.
        virtual uint32_t rand32() = 0;
        
        /// Reset the generator's state with the bits of `seed` as its source of entropy.
        virtual void rseed(uint64_t seed) = 0;
        
        /**
         * Produces a (pseudo-) random `T` with uniform distribution over the 
         * unit interval, if the unit interval is populated by `T`; or with uniform 
         * distribution over the entire space of possible `T` otherwise.
         * 
         * For example, `rand<float>()` returns a random float between 0 and
         * 1.0, while `rand<unsigned int>()` returns a random `uint` between 0 and `UINT_MAX`. 
         */
        template <typename T> T rand();
        
        /**
         * Produces a (pseudo-) random `T` with uniform distribution between
         * the `T` analogue of `0` and `max`.
         * 
         * @param max Upper bound of possible `T` samples.
         */
        template <typename T> T rand(T max);
        
        /**
         * Produces a (pseudo-) random `T` with uniform distribution between
         * `lo` and `hi`.
         * @param lo Lower bound for possible `T` samples.
         * @param hi Upper bound for possible `T` samples.
         */
        template <typename T> T rand(T lo, T hi);
        
    protected:
        uint32_t     _bitpool;
        unsigned int _bitsleft;
    };
}
#endif /* RANDOM_H_ */
