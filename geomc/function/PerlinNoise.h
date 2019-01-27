/*
 * PerlinNoise.h
 *
 *  Created on: Apr 9, 2011
 *      Author: tbabb
 */

#ifndef PERLINNOISE_H_
#define PERLINNOISE_H_

#include <boost/shared_array.hpp>
#include <boost/integer.hpp>
#include <limits>

#include <geomc/function/Utils.h>
#include <geomc/linalg/Vec.h>
#include <geomc/function/functiondetail/PerlinDetail.h>


#define PERLIN_NUM_GRADIENTS (0x100)
#define PERLIN_MODULO_MASK   (0xff)


namespace geom {

/**
 * @ingroup function
 * @brief Real-valued smooth noise over `N` dimensions.
 * 
 * Perlin noise has _O(2<sup>N</sup>)_ time cost.
 */
template <typename T, index_t N>
class PerlinNoise {
public:
    
    // define the grid type to be at least large enough to represent all
    // the integral values that can be represented by point_t.
    typedef PointType<index_t,N>        gridtype;
    typedef PointType<T,N>              pointtype;
    typedef typename gridtype::point_t  grid_t;
    typedef typename pointtype::point_t point_t;
    
    /**********************************
     * Structors                      *
     **********************************/
    
    /**
     * Construct a new perlin noise object with the default (non-reentrant) random number generator.
     */
    PerlinNoise():
            gradients(new point_t[PERLIN_NUM_GRADIENTS]),
            p(new index_t[PERLIN_NUM_GRADIENTS]) {
        detail::_ImplPerlinInit<T,N>::init(this, getRandom(), PERLIN_NUM_GRADIENTS);
    }
    
    /**
     * Construct a new perlin noise object with the supplied random number generator.
     * @param rng A source of random bits.
     */
    PerlinNoise(Random *rng):
            gradients(new point_t[PERLIN_NUM_GRADIENTS]),
            p(new index_t[PERLIN_NUM_GRADIENTS]) {
        detail::_ImplPerlinInit<T,N>::init(this, rng, PERLIN_NUM_GRADIENTS);
    }
    
    /**********************************
     * Expr Functions                 *
     **********************************/
    
    /**
     * Evaluate the noise at `pt`.
     */
    T eval(point_t pt) {
        point_t p0f = std::floor(pt);
        grid_t  p0  = (grid_t)p0f;    //grid pt
        point_t p_modulus = pt - p0f; //pos. within grid
        const index_t corners = 1 << N;
        T planeVals[corners];
        for (index_t c = 0; c < corners; c++) {
            grid_t cPt;
            for (index_t axis = 0; axis < N; axis++) {
                gridtype::iterator(cPt)[axis] = (c & (1 << axis)) != 0;
            }
            const point_t &grid_grad = get_grid_gradient(p0 + cPt);
            planeVals[c] = pointtype::dot((p_modulus - (point_t)cPt), grid_grad);
        }
        
        for (index_t axis = 0; axis < N; axis++) {
            T t = fade(p_modulus[axis]); //interpolation parameter
            for (index_t pair = 0; pair < (1 << (N-axis-1)); pair++) {
                planeVals[pair] = mix(t, planeVals[pair*2], planeVals[pair*2+1]);
            }
        }
        
        return planeVals[0];
    }
    
    /**********************************
     * Other Functions                *
     **********************************/
    
    /**
     * Evaluate the gradient of the noise function via discrete difference.
     * 
     * Has _2N_ cost over `eval()`.
     * @param pt
     * @return 
     */
    point_t gradient(point_t pt) {
        point_t g;
        //central difference method
        for (index_t axis = 0; axis < N; axis++) {
            point_t dx;
            pointtype::iterator(dx)[axis] = PERLIN_EPSILON;
            pointtype::iterator(g)[axis]  = eval(pt + dx) - eval(pt - dx);
        }
        return g;
    }
    
protected:
    
    inline const point_t& get_grid_gradient(const grid_t &pt) {
        index_t idx = 0;
        for (index_t axis = 0; axis < N; axis++) {
            idx = p[(idx + gridtype::iterator(pt)[axis]) & ((index_t)PERLIN_MODULO_MASK)];
        }
        return gradients[idx];
    }
    
    static inline T fade(T t) {
        return t*t*t*(t*(t*6 - 15) + 10);
    }
    
    void init(Random *rng) {
        detail::_ImplPerlinInit<T,N>::init(this, rng);
    }
    
    boost::shared_array<point_t> gradients;
    boost::shared_array<index_t> p;
    
    static T PERLIN_EPSILON;
    
    friend class detail::_ImplPerlinInit<T,N>;
};

//empirically tested for T = double (?)
//might be a good idea to specialize this for other T.
template <typename T, index_t N> T PerlinNoise<T,N>::PERLIN_EPSILON = 0.00001;

}

#endif /* PERLINNOISE_H_ */
