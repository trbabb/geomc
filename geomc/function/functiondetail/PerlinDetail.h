/*
 * PerlinDetail.h
 *
 *  Created on: Mar 23, 2013
 *      Author: tbabb
 */

#ifndef PERLINDETAIL_H_
#define PERLINDETAIL_H_

#include <geomc/random/Random.h>
#include <geomc/random/RandomTools.h>
#include <geomc/function/FunctionTypes.h>

namespace geom {
namespace detail {

// all of these bullshit helper classes are just c++'s utterly
// retarded, roundabout way of providing the functionality of a
// templated conditionally-compiled code block. If those were 
// just implemented directly, then this code wouldn't look so fucking 
// ugly. A pox on bjarne stroustrup for cursing humanity with his hideous 
// troll language.

template <typename T, int N>
class _ImplPerlinInit {
public:
    static void init(PerlinNoise<T,N> *pln, Random *rng, index_t n) {
        Sampler<T> rndVecs = Sampler<T>(rng);
        for (index_t i = 0; i < n; i++){
            pln->gradients[i] = rndVecs.template unit<N>();
            pln->p[i] = i;
        }
        permute<index_t>(pln->p.get(), n, *rng);
    }
};

template <typename T>
class _ImplPerlinInit<T,1> {
public:
    static void init(PerlinNoise<T,1> *pln, Random *rng, index_t n) {
        for (index_t i = 0; i < n; i++){
            pln->gradients[i] = rng->rand(-1,1); //that the right range? 
            pln->p[i] = i;
        }
        permute<index_t>(pln->p.get(), n, *rng);
    }
};

}; // namespace detail
}; // namespace geom

#endif /* PERLINDETAIL_H_ */
