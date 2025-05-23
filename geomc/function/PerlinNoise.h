#pragma once
/*
 * PerlinNoise.h
 *
 *  Created on: Apr 9, 2011
 *      Author: tbabb
 */

#include <utility>

#include <geomc/function/Utils.h>
#include <geomc/random/SampleVector.h>
#include <geomc/random/DenseDistribution.h>
#include <geomc/linalg/Vec.h>

namespace geom {

/**
 * @ingroup function
 * @brief Real-valued smooth noise over `N` dimensions.
 * 
 * Perlin noise has _O(2<sup>N</sup>)_ time cost.
 */
template <typename T, index_t N>
class PerlinNoise : public Dimensional<T,N> {
    static constexpr size_t N_GRADIENTS = 0x100;
    
public:
    
    using gridtype  = PointType<index_t,N>;
    using pointtype = PointType<T,N>;
    using grid_t    = typename gridtype::point_t;
    using typename Dimensional<T,N>::point_t;
    
private:
    std::shared_ptr<point_t[]> _gradients;
    
public:
    
    /**
     * Construct a new perlin noise object with a `std::random_device` as a source of
     * random bits.
     */
    PerlinNoise():PerlinNoise(std::random_device()) {}
    
    /**
     * Construct a new perlin noise object with the supplied random number generator.
     * @param rng A source of random bits.
     */
    template <typename Generator>
    PerlinNoise(Generator& rng): _gradients(new point_t[N_GRADIENTS]) {
        if constexpr (N == 1) {
            // gradients should be in the range -1,1
            DenseUniformDistribution<T> dist(-1, 1);
            for (index_t i = 0; i < N_GRADIENTS; i++) {
                _gradients[i] = dist(rng);
            }
        } else {
            // gradients are random unit vectors
            for (index_t i = 0; i < N_GRADIENTS; i++) {
                _gradients[i] = random_unit<T,N>(rng);
            }
        }
    }
    
    T operator()(point_t pt) const {
        return eval(pt);
    }
    
    /**
     * Evaluate the noise at `pt`.
     */
    T eval(point_t pt) const {
        const index_t corners = 1 << N;
        T f_x[corners];
        
        // we use floor() because a bare cast will round toward zero:
        point_t p0f   = std::floor(pt);
        grid_t  p0    = (grid_t)p0f;  // grid pt
        point_t p_mod = pt - p0f;     // pos. within grid
        
        // initialize our interpolants with (hyper-)planar values emanating from the corners.
        for (index_t c = 0; c < corners; c++) {
            grid_t cnr = corner(c);
            const point_t& grid_grad = get_grid_gradient(p0 + cnr);
            f_x[c] = pointtype::dot((p_mod - (point_t)cnr), grid_grad);
        }
        
        // we start by interpolating pairs of values along an arbitrary axis (axis 0).
        // those interpolated values will in turn be interpolated along the next axis,
        // as in bi-/trilinear interpolation. each interpolation reduces the number 
        // of values by half, until we have one value.
        
        for (index_t axis = 0; axis < N; axis++) {
            T x = pointtype::iterator(p_mod)[axis];
            T t = fade(x);
            for (index_t pair = 0; pair < (1 << (N - axis - 1)); pair++) {
                f_x[pair] = mix(t, f_x[pair * 2], f_x[pair * 2 + 1]);
            }
        }
        
        return f_x[0];
    }
    
    /**
     * Evaluate the gradient of the noise function at `pt`. 
     * 
     * @param pt Location at which to sample the noise function.
     * @return A pair of (`noise(x)`, `gradient(noise(x))`).
     */
    std::pair<T, point_t> gradient(point_t pt) const {
        const index_t corners = 1 << N;
        T           f[corners]; // values
        point_t df_dx[corners]; // gradients
        
        point_t p0f   = std::floor(pt);
        grid_t  p0    = (grid_t)p0f;
        point_t p_mod = pt - p0f;
        
        // useful facts:
        //      mix(t, a, b) = (1 - t) * a + t * b
        // d/dt mix(t, a, b) = b - a
        // d/da mix(t, a, b) = 1 - t
        // d/db mix(t, a, b) = t
        
        // the multivariate chain rule:
        // df/dt f(a(t), b(t), c(t)) = 
        //         df_da(a(t), ...)      * da_dt(t) + 
        //         df_db(..., b(t), ...) * db_dt(t) + ...
        
        // start with planar functions from each corner;
        // remember the derivatives of these functions:
        for (index_t c = 0; c < corners; ++c) {
            grid_t cnr = corner(c);
            df_dx[c]   = get_grid_gradient(p0 + cnr);
            f[c]       = pointtype::dot((p_mod - (point_t)cnr), df_dx[c]);
        }
        
        // now start interpolating these functions, updating the derivative along the way
        for (index_t axis = 0; axis < N; axis++) {
            T     x = pointtype::iterator(p_mod)[axis];
            T     t = fade(x);
            T dt_dx = dfade_dt(x); // (1) note that this is zero along other axes,
                                   // as this is only a fuction of x[axis].
            
            // as in `eval()`, we begin by interpolating pairs of points along the x axis
            // and reducing the number of values by half with successive interpolations.
            for (index_t pair = 0; pair < (1 << (N - axis - 1)); pair++) {
                index_t k = pair * 2;
                T a = f[k];
                T b = f[k + 1];
                point_t da_dx = df_dx[k];
                point_t db_dx = df_dx[k + 1];
                
                // the expression to be differentiated:
                f[pair] = mix(t, a, b);
                
                // apply the multivariate chain rule:
                point_t w;
                // dt/dx is only nonzero along the current axis:
                T* w_i  = pointtype::iterator(w) + axis;
                //        dmix_dt(t, ...) * dt/dx
                *w_i   += (b - a) * dt_dx;
                //        dmix_da(..., a, ...) * da/dx
                w      += (1 - t) * da_dx;
                //        dmix_db(..., ..., b) * db/dx
                w      += t * db_dx;
                df_dx[pair] = w;
            }
        }
        
        return std::pair<T, point_t>(f[0], df_dx[0]);
    }
    
protected:
    
    const point_t& get_grid_gradient(const grid_t& pt) const {
        uint64_t idx = 0;
        auto m = gridtype::iterator(pt);
        for (index_t i = 0; i < N; ++i) {
            // do an iterated linear congruential scramble, using Knuth's constants:
            uint64_t k = static_cast<uint64_t>(m[i]);
            idx = 6364136223846793005ULL * (k + idx) + 1442695040888963407ULL;
        }
        idx = static_cast<index_t>(positive_mod<uint64_t>(idx, N_GRADIENTS));
        return _gradients[idx];
    }
    
    static inline grid_t corner(index_t i) {
        grid_t cnr;
        for (index_t axis = 0; axis < N; axis++) {
            gridtype::iterator(cnr)[axis] = (i & (1 << axis)) != 0;
        }
        return cnr;
    }
    
    static inline T fade(T t) {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }
    
    static inline T dfade_dt(T t) {
        T w = t - 1;
        return 30 * w * w * t * t;
    }
};

}
