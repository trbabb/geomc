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

/*
 * random has some design problems:
 *   1) Templated types (Complex<T>, e.g.) cannot have rand() specializations,
 *      because c++ (for no good reason) does not allow partially-specialized 
 *      template functions.
 *   2) There wouldn't be a way to make a drop-in Random class which, e.g., implements
 *      a low discrepancy sequence, because rand() isn't virtual.
 *   3) Picking from distributions is mostly unaddressed.
 *      Note that Sampler is, in a way, a tool for picking from distributions over
 *      vectors, though it might also be a sensible thing to pass other 
 *      distributions to the sampler for radial density, angular density, etc. etc.
 
 *      I should probably lookup std::rand and see if maybe that's actually better.
 *      > note that std::rand is c++11 only.
 * 
 * requirements for new design: 
  - picking from uniform distributions is easy.
  - adding new types of objects to pick is easy.
  - it is possible to partially specialize picking of templated objects.
  - objects can be picked from a source of entropy in a predictable sequence.
  - low discrepancy sequences have a natural representation.
  ~ would be nice if swapping out other distributions for uniforms were 
    possible.
 * 
 * proposal:
 *  rand<type,distribution=UNIT_INTERVAL>() calls RandomHelper<type>::rand()
 *  rand<type>(hi)    // half-interval
 *  rand<type>(lo,hi) // full-interval
 * 
 * (would functors help us?)
 * 
 * other distributions such as: normal<params> or statified<params> or hammersley<dimensions>
 * So what you need is:
 *   a way to slip something in between the unit-interval picker and the final result
 *   which:
 *      - doesn't make basic picking slower
 *      - doesn't make uniform interval choices slower
 *      - doesn't make unit and interval choices more difficult or cumbersome to use.
 * 
 * There is something to be said for simply transforming floats-- assumed to 
 * be uniform-- to shapes with simple functions that are agnostic of the generator.
 * This is certainly better for stratified/LD sampling.
 * 
 * The one thing I forbid is constructing a new object specifically to pick from a 
 * uniform distribution. That's too verbose.
 * 
 * what we have is
 *   entropy -> PRNG -> (random bits) -> uniform real^n -> distribution
 *    kb?       MT?                      statified?        uniform?
 *    time?     LC?                      Downey?           normal?
 *    network?  yarrow?                  Hammersley?       poisson?
 *    device?   etc.                     random?           etc.
 *    hard seed?
 *    etc.
 * 
 * where each of those parts could be interchanged.
 * 
 * entropy:
 *    PRNG has a method, seed(), which takes a char* and a count, for entropy
 *    input. maybe this is different for each PRNG (e.g. for Yarrow, which
 *    requires an estimate of entropy bits. at the very least,
 *    each PRNG should take a single 32 or 64 bit seed. beyond that, modularity
 *    is not important.
 *       note that Hammersley (?) is completely deterministic and doesn't take a 
 *       seed.
 * PRNG:
 *    this is pretty straightforward, a PRNG's only job is to output random bits.
 *    so having a method to fill an array with bits is maybe a useful thing,
 *    rather than mandating 32 or 64 bits at a time. But providing these methods
 *    shall be provided for bare minimum interoperability. the rand<32> and rand(ct)
 *    methods may be implemented one-in-terms-of-the-other in whichever direction
 *    is convenient for the implementation.
 * Uniform real:
 *    perhaps this is a template param of rand<T>()? a class to transform random
 *    bits to (0,1)^n. 
 *    - may need more than a class name-- instances might be necessary.
 *      e.g. a stratified sampler with parameterized strata counts.
 *    - if not bundled with PRNG, hard to make modular/swappable.
 *    - if bundled with PRNG, stateful. 
 *    - Some algorithms are not parameterized (Downey) and some are (stratified).
 *    - some care about how many samples are to be drawn, and some don't.
 *      This makes encapsulation harder.
 *    This is the object that we should pass around and call Random.
 *    Each Random is templated over a single type of object.
 *    Distributions over the reals^N are handled by a specialization of 
 *    Random over vectors. 
 *    > i.e., we make a Hammersley<T,N> which is an instance of Random<Vec<T,N>>?
 *    - yes, I think this is right. but the last question is the virtual indirection 
 *      we'll be adding. can we remove that without adding too much complexity?
 * Distribution:
 *    Drawing from the uniform distribution should not be hard or require 
 *    constructing objects. This should probably be its own class,
 *    and does not need to be "drop-in", in that the client of the random
 *    numbers should be in control of the distribution, and will likely need to 
 *    know if it is changing. if modularity is desired, the client should explicitly
 *    provide/expose it.
 *    
 *    basically, code is likely to make the assumption that the chosen numbers are
 *    from a uniform distribution, and we can't allow that to be easily broken by
 *    passing a magically-wrapped Random with some j. random distribution swapped in.
 * 
 *    As such we don't need to implement distributions just yet; they're pretty 
 *    orthogonal to the rest of the number generation problem.
 * 
 * note that really, distributions (gaussian, e.g.) may not pick their numbers 
 * by transforming a uniform real. they might pick based on the pseudorandom 
 * bits directly. This is especially true if we model LD sequences as distributions.
 * 
 * There is nothing prohibiting us from nesting distributions, e.g. passing a 
 * uniform real distribution to a gaussian distribution. So maybe we should
 * make 'uniform real' a special type of distribution which is virtual, and that's 
 * what we normally pass around.
 */

/**
 * @defgroup random
 * @brief Random number and object generation.
 * 
 * To obtain random numbers, choose and construct a random number generator,
 * then call `rand()` with the desired type as a template parameter:
 * 
 *     Random *rng = new MTRand(mySeed);
 *     float f = rng->rand<float>(); // a random number between 0 and 1.0 
 * 
 * Or use `geom::getRandom()` for a quick-and-dirty (not reentrant or threadsafe!) 
 * generator instance:
 * 
 *     float f = getRandom()->rand<float>();
 * 
 * A note about random ranges
 * --------------------------
 * 
 * For integer random numbers in a range, it is far preferable to use the provided
 * methods, rather than the "standard" `rand() % rangeMax` strategy, which is not 
 * a uniform distribution, most severely when `rangeMax` is a large fraction of 
 * `RAND_MAX`. 
 * 
 * The provided methods of generating floating point numbers are also generally 
 * superior to the (again) "standard" method of calculating `rand() / RAND_MAX`,
 * which excludes about 93% of representable numbers between 0 and 1, particularly
 * in the regions near 0. `geom::Random::rand(float)` addresses this by carefully
 * choosing the bits of the generated number explicitly. The method is fast, makes
 * efficient use of the supplied entropy, and produces a uniform distribution. (Note
 * that C++11's `<random>` library does *not* address this issue).
 * 
 */

// Notice: This #define is referred to in other places, namely Dual.h.
//         This is part of a mechanism to bring in code which is common
//         to both `random` and `function` iff both libraries are in use.

#ifndef GEOMC_RANDOM_H_
#define GEOMC_RANDOM_H_

#include <geomc/geomc_defs.h>
#include <geomc/random/randomdetail/RandomImpl.h>


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
        
        // todo: make the base template of rand<>() call randomimpl<T>
        // user specializations must go in randomimpl.
        // add specialization declarations to the random class here, in the
        // header. this will prevent the base template from masking the
        // specializations, and might also prevent doxygen from getting confused.
        template <typename T> T rand() {
            return detail::RandomImpl<T>::rand(this);
        }
        
        /**
         * Produces a (pseudo-) random `T` with uniform distribution between
         * the `T` analogue of `0` and `max`.
         * 
         * @param max Upper bound of possible `T` samples.
         */
        template <typename T> T rand(T max) {
            return detail::RandomImpl<T>::rand(this, max);
        }
        
        /**
         * Produces a (pseudo-) random `T` with uniform distribution between
         * `lo` and `hi`.
         * @param lo Lower bound for possible `T` samples.
         * @param hi Upper bound for possible `T` samples.
         */
        template <typename T> T rand(T lo, T hi) {
            return detail::RandomImpl<T>::rand(this, lo, hi);
        }

        
    protected:
        uint32_t     _bitpool;
        unsigned int _bitsleft;
    };
}


// oddly enough, this goes at the bottom, because the inner template function
// needs to be visible before it can be fully specialized. god, c++ is dumb

#include <geomc/random/randomdetail/RandomSpecialize.h>

// same here. RandomDual needs to see the rand<>() specialization defs so it can
// call them. So fragile. Ick.

#ifdef GEOMC_DUAL_H
#include <geomc/function/functiondetail/RandomDual.h>
#endif


#endif /* GEOMC_RANDOM_H_ */
