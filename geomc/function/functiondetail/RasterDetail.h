/*
 * RasterDetail.h
 *
 *  Created on: Mar 16, 2013
 *      Author: tbabb
 */

#ifndef RASTERDETAIL_H_
#define RASTERDETAIL_H_

#include <geomc/function/FunctionTypes.h>

namespace geom {
namespace detail {

/*************************
 * Edge behaviors        *
 *************************/

template <EdgeBehavior Edge>
class _ImplEdge {
    // no contents.
};

template <>
class _ImplEdge<EDGE_CLAMP> {
public:
    static inline index_t coord(index_t c, index_t max) {
        return std::max(std::min(c,max-1),(index_t)0);
    }
};

template <>
class _ImplEdge<EDGE_PERIODIC> {
public:
    static inline index_t coord(index_t c, index_t max) {
        return positive_mod(c,max);
    }
};

template <>
class _ImplEdge<EDGE_MIRROR> {
public:
    static inline index_t coord(index_t c, index_t max) {
        c = std::abs(c);
        max -= 1;
        int cyc = c/max;
        int i   = c%max;
        return ((cyc & 1) == 0)?(i):(max-i);
    }
};

template <>
class _ImplEdge<EDGE_CONSTANT> {
public:
    static inline index_t coord(index_t c, index_t max) {
        return c;
    }
};

// the following is completely unreadable
// because c++ is utterly moronic.
// this should all be in the main class, but
// member functions cannot be specialized in-place
// because c++ is, again, too damn stupid to
// perform the simple substitution transformtaion
// applied here. therefore we make a fuckton of 
// opaque helper classes, one for pretty much every 
// single member function. fuck me. 
// fuck c++. fuck bjarne stroustrup.

/*************************
 * Sampling behaviors    *
 *************************/


template <typename I, typename O, index_t N, index_t Channels, EdgeBehavior Edge, Interpolation Interp>
class _ImplRasterSample {
    // pass
};

template <typename I, typename O, index_t N, index_t Channels, EdgeBehavior Edge>
class _ImplRasterSample<I,O,N,Channels,Edge,INTERP_NEAREST> {
public:
    typedef typename Raster<I,O,N,Channels>::coord_t  coord_t;
    typedef typename Raster<I,O,N,Channels>::grid_t   grid_t;
    typedef typename Raster<I,O,N,Channels>::sample_t sample_t;
    
    static inline sample_t sample(const Raster<I,O,N,Channels> *r, const coord_t &pt) {
        grid_t gridPt = grid_t(pt + coord_t(0.5));
        return r->template sample_discrete<Edge>(gridPt);
    }
};

template <typename I, typename O, index_t N, index_t Channels, EdgeBehavior Edge>
class _ImplRasterSample<I,O,N,Channels,Edge,INTERP_LINEAR> {
public:
    typedef typename Raster<I,O,N,Channels>::coord_t  coord_t;
    typedef typename Raster<I,O,N,Channels>::grid_t   grid_t;
    typedef typename Raster<I,O,N,Channels>::sample_t sample_t;
    
    static inline sample_t sample(const Raster<I,O,N,Channels> *r, const coord_t &pt) {
        sample_t buf[1<<N];
        grid_t gridPt = (grid_t)pt;
        coord_t s = pt - ((coord_t)gridPt);
        
        // copy surrounding 2^N sample pts into a contiguous buffer
        r->template copy<Edge>(buf, Rect<int,N>(gridPt, gridPt + grid_t(2)));
        
        return interp_linear(PointType<I,N>::iterator(s), buf, N);
    }
};

template <typename I, typename O, index_t N, index_t Channels, EdgeBehavior Edge>
class _ImplRasterSample<I,O,N,Channels,Edge,INTERP_CUBIC> {
public:
    typedef typename Raster<I,O,N,Channels>::coord_t  coord_t;
    typedef typename Raster<I,O,N,Channels>::grid_t   grid_t;
    typedef typename Raster<I,O,N,Channels>::sample_t sample_t;
    
    static inline sample_t sample(const Raster<I,O,N,Channels> *r, const coord_t &pt) {
        sample_t buf[1<<(2*N)];
        grid_t gridPt = (grid_t)pt;
        coord_t s = pt - ((coord_t)gridPt);
        
        // copy surrounding 4^N sample pts into a contiguous buffer
        r->template copy<Edge>(buf, Rect<int,N>(gridPt - grid_t(1), gridPt + grid_t(3)));
        
        return interp_cubic(PointType<I,N>::iterator(s), buf, N);
    }
};

/*************************
 * Indexing strategy     *
 *************************/

template <typename grid_t, index_t N, EdgeBehavior Edge>
class _ImplRasterIndex {
public:
    
    static inline index_t index(const grid_t &extent, const grid_t &c) {
        index_t dim = 1;
        index_t idx = 0;
        for (int i = 0; i < N; i++){
            index_t x = detail::_ImplEdge<Edge>::coord(c[i], extent[i]);
            idx += x*dim;
            dim *= extent[i];
        }
        return idx;
    }
};

template <typename grid_t, EdgeBehavior Edge>
class _ImplRasterIndex<grid_t, 1, Edge> {
public:
    
    static inline index_t index(const grid_t &extent, const grid_t &c) {
        return detail::_ImplEdge<Edge>::coord(c, extent);
    }
};

/*************************
 * Volume calc           *
 *************************/

template <index_t N>
index_t array_product(const index_t *start) {
    index_t prod = 1;
    const index_t *end = start+N;
    for (; start < end; start++){
        prod *= *start;
    }
    return prod;
}

template <>
inline index_t array_product<1>(const index_t *start) {
    return *start;
}


}; // namespace detail
}; // namespace geom

#endif /* RASTERDETAIL_H_ */
