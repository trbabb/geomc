/*
 * Raster.h
 *
 *  Created on: Mar 12, 2013
 *      Author: tbabb
 */

#ifndef RASTER_H_
#define RASTER_H_

#include <algorithm>
#include <boost/shared_array.hpp>

#include <geomc/linalg/Vec.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/GridIterator.h>
#include <geomc/function/Utils.h>
#include <geomc/function/functiondetail/RasterDetail.h>

namespace geom {

//TODO: remove cluttering templated sample functions? (15% speed difference, though)
//TODO: "lobe-yness" of interp is subject to parameterization
//TODO: really, interpolation and edge behavior should be objects.
//      e.g. abyss color belongs to abyss edge behavior only.
//      - could parameterize in the Raster object. i.e. r.setEdgeBehavior()
//        or in the call site, i.e. r.sample(pt, EdgeBehavior(...));
//TODO: make work for dynamic.
//      (specialize entire class?)
//      or, class RasterHandle<I,O,M,N> : public virtual AnyRaster
//TODO: add:
//      - resample<MAG_FILTER,MIN_FILTER>(grid)
//      - resample<...>(affineTransform)
//      - resample<...>(raster<grid_t->in_t> pts)
//      see http://entropymine.com/imageworsener/resample/


/*************************
 * Raster class          *
 *************************/

// should be re-templated:
// - input dimensions?
// - point type
// - point dimension
// - out type (defaults to point type)
//   > can't template over input type at call site because compiler will get confused by:
//       template <typename T,K>
//       void sample(const typename PointType<T,K>::point_t &pt)

/**
 * @ingroup function
 * @brief An M-dimensional grid of interpolated data which can be continuously sampled.
 * 
 * @tparam I Domain data type.
 * @tparam O Range data type.
 * @tparam M Domain dimension.
 * @tparam N Range dimension.
 * 
 * In other words, a function of I<sup>M</sup> &rarr; O<sup>N</sup>.
 * 
 * Some functions are overloaded to have some arguments as template parameters.
 * This allows some low-level sampling code to be inlined which can improve
 * speed by as much as 15%. The traditional function-argument versions are available
 * for when the parameters must be chosen at runtime.
 */
template <typename I, typename O, index_t M, index_t N>
class Raster {
public:
    
    /// Type of sample location. `I` if `M` is 1, otherwise `Vec<I,M>`.
    typedef typename PointType<I,M>::point_t        coord_t;
    /// Type for indexing data, i.e. a grid coordinate. `index_t` if `M` is 1, `Vec<index_t,M>` otherwise.
    typedef typename PointType<index_t,M>::point_t  grid_t;
    /// Type of resultant data. `O` if `N` is 1, otherwise `Vec<O,N>`.
    typedef typename PointType<O,N>::point_t sample_t;
    
protected:
    grid_t  extent;
    index_t size;
    boost::shared_array<O> data;
    sample_t abyss;
    
public:
    
    ////////// Structors //////////
    
    /// Construct a Raster with `dims` elements along each axis.
    Raster(const grid_t &dims):
            extent(dims),
            size(detail::array_product<M>(PointType<index_t,M>::iterator(dims)) * N),
            data(new O[size]){
        std::fill(data.get(), data.get()+size, 0);
    }
    
    ////////// Methods //////////
    
    /**
     * Set the sample value to be used for "out-of-bounds" samples when 
     * `EDGE_CONSTANT` behavior is used.
     */
    void setAbyss(sample_t val) {
        abyss = val;
    }
    
    /**
     * Set the value of the raster at the point `pt`. 
     * @param pt Coordinate of datapoint to set.
     * @param val New value of datapoint.
     */
    void set(const grid_t &pt, const sample_t &val) {
        if (contains(pt)){
            O *p = data.get() + this->template index<EDGE_CONSTANT>(pt);
            //xxx no worky with dynamic
            *(sample_t*)p = val;
        }
    }
    
    /**
     * Copy the data in `region` into the array `dest`, in row-major 
     * (first dimension consecutive) order. Templated over edge sampling
     * behavior.
     * 
     * @tparam Edge Edge sampling behavior.
     * @param dest Destination array.
     * @param region Region of this Raster to copy from.
     */
    template <EdgeBehavior Edge>
    void copy(sample_t *dest, const Rect<index_t,M> &region) const {
        GridIterator<index_t,M> i = GridIterator<index_t,M>(region);
        GridIterator<index_t,M> end = i.end();
        for (; i != end; ++i, ++dest){
            *dest = this->template sample_discrete<Edge>(*i);
        }
    }
    
    /**
     * Copy the data in `region` into the array `dest`, in row-major
     * (first dimension consecutive) order. Dynamic edge sampling behavior.
     * 
     * @param dest Destination array.
     * @param region Region of this raster to copy from.
     * @param edge Edge sampling behavior.
     */
    void copy(sample_t *dest, const Rect<index_t,M> &region, EdgeBehavior edge) const {
        GridIterator<index_t,M> i = GridIterator<index_t,M>(region);
        GridIterator<index_t,M> end = i.end();
        for (; i != end; ++i, ++dest){
            *dest = this->sample_discrete(*i, edge);
        }
    }
    
    /**
     * Sample this Raster at `pt`.
     * 
     * @tparam Edge Edge sampling behavior.
     * @tparam Interp Sample interpolation strategy.
     * @param pt Point to sample.
     * @return Sampled data.
     */
    template <EdgeBehavior Edge, Interpolation Interp>
    inline sample_t sample(const coord_t &pt) const {
        return detail::_ImplRasterSample<I,O,M,N,Edge,Interp>::sample(this, pt);
    }
    
    /**
     * Sample this Raster at `pt`.
     * 
     * @param pt Point to sample.
     * @param edge Edge sampling behavior.
     * @param interp Sample interpolation strategy.
     * @return Sampled data.
     */
    sample_t sample(const coord_t &pt, EdgeBehavior edge=EDGE_CLAMP, Interpolation interp=INTERP_LINEAR) const {
        grid_t gridPt;
        
        switch (interp){
        case INTERP_NEAREST:
            gridPt = grid_t(pt + coord_t(0.5));
            return sample_discrete(gridPt, edge);
        case INTERP_LINEAR:
        {
            sample_t buf_l[1<<M];
            gridPt = (grid_t)pt;
            coord_t s = pt - ((coord_t)gridPt);
            
            // copy surrounding 2^M sample pts into a contiguous buffer
            copy(buf_l, Rect<index_t,M>(gridPt, gridPt + grid_t(2)), edge);
            
            return interp_linear(s.begin(), buf_l, M);
        }
        case INTERP_CUBIC:
        {
            sample_t buf_c[1<<(2*M)];
            gridPt = (grid_t)pt;
            coord_t s = pt - ((coord_t)gridPt);
            
            // copy surrounding 4^M sample pts into a contiguous buffer
            copy(buf_c, Rect<index_t,M>(gridPt - grid_t(1), gridPt + grid_t(3)), edge);
            
            return interp_cubic(s.begin(), buf_c, M);
        }
        default:
            return abyss;
        }
    }
    
    /**
     * Retrieve the data at discrete grid location `pt`. 
     * 
     * @tparam Edge sampling behavior.
     * @param pt Grid location to sample.
     * @return Data at `pt`.
     */
    template <EdgeBehavior Edge>
    sample_t sample_discrete(const grid_t &pt) const {
        index_t offs;
        switch (Edge){
            case EDGE_CONSTANT:
                if (contains(pt)){
                    offs = this->template index<EDGE_CONSTANT>(pt);
                } else {
                    return abyss;
                }
                break;
            default:
                offs = this->template index<Edge>(pt);
                break;
        }
        return PointType<O,N>::from_ptr(data.get() + offs);
    }
    
    /**
     * Retrieve the data at discrete grid location `pt`.
     * 
     * @param pt Grid location to sample.
     * @param edge Edge sampling behavior.
     * @return Data at `pt`.
     */
    sample_t sample_discrete(const grid_t &pt, EdgeBehavior edge) const {
        index_t offs;
        
        switch (edge){
        case EDGE_CONSTANT:
            if (!contains(pt)) return abyss;
            offs = this->template index<EDGE_CONSTANT>(pt);
            break;
        case EDGE_CLAMP:
            offs = this->template index<EDGE_CLAMP>(pt);
            break;
        case EDGE_PERIODIC:
            offs = this->template index<EDGE_PERIODIC>(pt);
            break;
        case EDGE_MIRROR:
            offs = this->template index<EDGE_MIRROR>(pt);
            break;
        default:
            return abyss;
        }
        
        return sample_t(data.get() + offs);
    }
    
    /**
     * @return `true` if the grid location `pt` is within the bounds of this
     * raster and has data, `false` otherwise.
     */
    bool contains(const grid_t &pt) const {
        typedef PointType<index_t,M> grid_info;
        
        const index_t *p_i = grid_info::iterator(pt);
        const index_t *e_i = grid_info::iterator(extent);
        for (index_t i = 0; i < M; i++){
            index_t c = p_i[i];
            if (c < 0 or c >= e_i[i]) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Size of this raster along each axis.
     */
    inline grid_t extents() const {
        return extent;
    }
    
    /**
     * Number of input dimensions.
     */
    inline index_t inputDimension() const {
        return M;
    }
    
    /**
     * Number of output dimensions.
     */
    inline index_t outputDimension() const {
        return N;
    }
    
protected:
    template <EdgeBehavior Edge>
    inline index_t index(const grid_t &c) const {
        // specialized where M=1
        return detail::_ImplRasterIndex<grid_t,M,Edge>::index(this->extent, c) * outputDimension();
    }
};

}; // namespace geom

#endif /* RASTER_H_ */
