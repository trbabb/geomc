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

#include <geomc/Utils.h>
#include <geomc/linalg/Vec.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/GridIterator.h>
#include <geomc/function/functiondetail/RasterDetail.h>

namespace geom {

//TODO: remove cluttering templated sample functions? (15% speed difference, though)
//TODO: "lobe-yness" of interp is subject to parameterization
//TODO: make work for dynamic.
//      (specialize entire class?)
//      or, class RasterHandle<I,O,N,Channels> : public virtual AnyRaster
//TODO: add:
//      - resample<MAG_FILTER,MIN_FILTER>(grid)
//      - resample<...>(affineTransform)
//      - resample<...>(raster<grid_t->in_t> pts)
//      see http://entropymine.com/imageworsener/resample/


/*************************
 * Raster class          *
 *************************/

template <typename I, typename O, index_t N, index_t Channels>
class Raster {
public:
    
    typedef typename PointType<I,N>::point_t        coord_t;
    typedef typename PointType<index_t,N>::point_t  grid_t;
    typedef typename PointType<O,Channels>::point_t sample_t;
    
protected:
    grid_t  extent;
    index_t size;
    boost::shared_array<O> data;
    sample_t abyss; //xxx no worky with dynamic
    
public:
    
    ////////// Structors //////////
    
    Raster(const grid_t &dims):
            extent(dims),
            size(detail::array_product<N>(PointType<index_t,N>::iterator(dims)) * Channels),
            data(new O[size]){
        std::fill(data.get(), data.get()+size, 0);
    }
    
    ////////// Methods //////////
    
    void setAbyss(sample_t val) {
        abyss = val;
    }
    
    void set(const grid_t &pt, const sample_t &val) {
        if (contains(pt)){
            O *p = data.get() + this->template index<EDGE_CONSTANT>(pt);
            //xxx no worky with dynamic
            *(sample_t*)p = val;
        }
    }
    
    //xxx not work with dynamic
    template <EdgeBehavior Edge>
    void copy(sample_t *dest, const Rect<index_t,N> &region) const {
        GridIterator<index_t,N> i = GridIterator<index_t,N>(region);
        GridIterator<index_t,N> end = i.end();
        for (; i != end; ++i, ++dest){
            *dest = this->template sample_discrete<Edge>(*i);
        }
    }
    
    void copy(sample_t *dest, const Rect<index_t,N> &region, EdgeBehavior edge) const {
        GridIterator<index_t,N> i = GridIterator<index_t,N>(region);
        GridIterator<index_t,N> end = i.end();
        for (; i != end; ++i, ++dest){
            *dest = this->sample_discrete(*i, edge);
        }
    }
    
    template <EdgeBehavior Edge, Interpolation Interp>
    inline sample_t sample(const coord_t &pt) const {
        return detail::_ImplRasterSample<I,O,N,Channels,Edge,Interp>::sample(this, pt);
    }
    
    sample_t sample(const coord_t &pt, EdgeBehavior edge=EDGE_CLAMP, Interpolation interp=INTERP_LINEAR) const {
        grid_t gridPt;
        
        switch (interp){
        case INTERP_NEAREST:
            gridPt = grid_t(pt + coord_t(0.5));
            return sample_discrete(gridPt, edge);
        case INTERP_LINEAR:
        {
            sample_t buf_l[1<<N];
            gridPt = (grid_t)pt;
            coord_t s = pt - ((coord_t)gridPt);
            
            // copy surrounding 2^N sample pts into a contiguous buffer
            copy(buf_l, Rect<index_t,N>(gridPt, gridPt + grid_t(2)), edge);
            
            return interp_linear(s.begin(), buf_l, N);
        }
        case INTERP_CUBIC:
        {
            sample_t buf_c[1<<(2*N)];
            gridPt = (grid_t)pt;
            coord_t s = pt - ((coord_t)gridPt);
            
            // copy surrounding 4^N sample pts into a contiguous buffer
            copy(buf_c, Rect<index_t,N>(gridPt - grid_t(1), gridPt + grid_t(3)), edge);
            
            return interp_cubic(s.begin(), buf_c, N);
        }
        default:
            return abyss;
        }
    }
    
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
        return PointType<O,Channels>::from_ptr(data.get() + offs);
    }
    
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
    
    bool contains(const grid_t &pt) const {
        typedef PointType<index_t,N> nfo;
        
        for (index_t i = 0; i < N; i++){
            index_t c = nfo::iterator(pt)[0];
            if (c < 0 or c >= nfo::iterator(extent)[i]) {
                return false;
            }
        }
        return true;
    }
    
    inline grid_t extents() const {
        return extent;
    }
    
    inline index_t numChannels() const {
        return Channels;
    }
    
protected:
    template <EdgeBehavior Edge>
    inline index_t index(const grid_t &c) const {
        // specialized where N=1
        return detail::_ImplRasterIndex<grid_t,N,Edge>::index(this->extent, c) * numChannels();
    }
};

}; // namespace geom

#endif /* RASTER_H_ */
