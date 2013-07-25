/*
 * Random.h
 * 
 * Abstract class for random number generation.
 * Implementations provide random bits by overriding rand32(), which should return 
 * no fewer than 32 random bits.
 * 
 * Client-friendly random numbers are exposed via Random::rand<T>(). Specializations
 * are provided for a selection of primitive T. Specializations for other T may be written,
 * without the need to create new subclasses.
 * 
 * 
 * The general contract of rand() which should be followed by all specializations:
 * 
 *   rand(void)
 *     - if T is discrete, than this function should return values evenly over the entire 
 *       space of possible T.
 *     - if T is a continuous type, then this function should return values evenly
 *       over an analogue of the interval [0,1], inclusive.
 * 
 *   rand(T max)
 *     - should return random values evenly distributed over all possible states between
 *       0 (or the T analogue of 0) and the value of <max>
 * 
 *   rand(T lo, T hi)
 *     - should return random values evenly distributed over all possible states 
 *       bounded by <lo> and <hi>
 * 
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

#ifndef RANDOM_H_
#define RANDOM_H_

#include <stdint.h>

namespace geom {
    class Random {
    public:

        Random(void);
        Random(uint64_t seed); //allow us to supply up to 64 bits of entropy
        virtual ~Random();

        //returns 32 random bits.
        virtual uint32_t rand32() = 0;
        virtual void rseed(uint64_t seed) = 0;
        
        template <typename T> T rand();
        template <typename T> T rand(T max);
        template <typename T> T rand(T lo, T hi);
        
    protected:
        uint32_t     _bitpool;
        unsigned int _bitsleft;
    };
}
#endif /* RANDOM_H_ */
