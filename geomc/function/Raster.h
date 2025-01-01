/*
 * Raster.h
 *
 *  Created on: Mar 12, 2013
 *      Author: tbabb
 */

#ifndef RASTER_H_
#define RASTER_H_

#include <algorithm>

#include <geomc/linalg/Vec.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/GridIterator.h>
#include <geomc/function/Utils.h>
#include <geomc/function/functiondetail/RasterDetail.h>

namespace geom {

//TODO: I question the design of this class.
//      all this edge behavior bullshit should be handled simply by creating
//      different instances with different behavior. They can share their memory
//      underneath; constructing one from the other should be cheap/easy.
//      > Use CRTP to preserve inline-ability?
//      > functions like integrate and resample will also be able to 
//        bypass the need for interp/edge params because it's internal state.
//      > but... need to pass around rasters without these stupid template params.
    
/*
 * Need two classes: Raster and Grid (or some clearer naming thereof). 
 * Grid always has array indexes and is not interpolated.
 * Grids support basic arithmetic ops.
 * Raster wraps a Grid, interpolates it, and defines a domain.
 * Should grids be agnostic of input/output dim? I.e. always functions of
 * T^n -> T? and simply provide functions for slicing differently?
 *   what was the reason we got rid of Raster<Vec<T,N>> ? did we actually try that?
 */
/*
 * new class design
 *   - base class which handles base data, extents, domain, and the exposure thereof
 *   - derived needs to be either static or dynamic for each of (edge, interp)
 *   - edge behavior (like abyss) invites class variants with possible members
 *   - edge behaviors know something about integration
 *   - interp_nearest has phase parameters (samples are centered or what?)
 */
//TODO: have a raster which can wrap an array owned by someone else.
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
//      - integrate(coord_t min, coord_t max)
//      - arithmetic operators.
//        x domain issues. what if domain mismatches?
//          could make an absurdly huge new domain if using 'union' behavior
//      see http://entropymine.com/imageworsener/resample/
//TODO: reshaping dimensions should be easy.


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
 * Coordinates
 * ===========
 * 
 * Rasters have two coordinate systems: Indexing space and sampling space. 
 * In indexing space, the stored sample values may be accessed and changed 
 * using integer grid coordinates, analogous to a multi-dimensional array. 
 * 
 * Rasters may also be sampled continuously, and intermediate values between sample
 * points will be generated via interpolation. The domain of sampling space may be 
 * specified by the user, and sampling coordinates will be renormalized to this region.
 * If no domain is specified, sampling space is equal to indexing space (with the
 * addition that intermediate values may be sampled). The user-specified domain
 * has no effect on indexing space.
 * 
 * Sampling functions
 * ================== 
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
    typedef typename PointType<O,N>::point_t        sample_t;
    
protected:
    grid_t    m_extent;
    index_t   m_size;
    std::shared_ptr<O[]> m_data;
    sample_t  m_abyss;
    Rect<I,M> m_domain;
    
public:
    
    ////////// Structors //////////
    
    /**
     * Construct a raster with `dims` elements along each axis. The sampling 
     * domain will be constructed such that the lower extreme is at the origin and
     * sample points are placed at integer coordinates.
     * 
     * @param dims Number of samples along each axis.
     */
    Raster(const grid_t &dims):
            m_extent(dims),
            m_size(detail::array_product<M>(PointType<index_t,M>::iterator(dims)) * N),
            m_data(new O[m_size]),
            m_domain((I)0, (coord_t)dims - (coord_t)1) {
        std::fill(m_data.get(), m_data.get()+m_size, 0);
    }
    
    /**
     * Construct a raster with `dims` elements along each axis.
     * 
     * @param dims Number of samples along each axis.
     * @param domain Desired boundary of the data region. The lower and upper 
     * coordinates correspond to the exact coordinates of the most extreme data 
     * points along each axis.
     */
    Raster(const grid_t &dims, const Rect<I,M> &domain):
            m_extent(dims),
            m_size(detail::array_product<M>(PointType<index_t,M>::iterator(dims)) * N),
            m_data(new O[m_size]),
            m_domain(domain) {
        std::fill(m_data.get(), m_data.get()+m_size, 0);
    }
    
    /** 
     * Construct a Raster with `dims` elements along each axis.
     * Fill the raster with the data in `src_data`. `src_data` must have enough 
     * elements to fill the raster completely; in other words, the length of `src_data`
     * must be the product of the elements of `dims`. The domain will be constructed
     * to align with the sample grid coordinates; in other words, the lowest 
     * extreme is at the origin and sample points are placed at integer coordinates.
     * 
     * @param dims Number of samples along each axis.
     * @param src_data Data to copy into this raster.
     */
    Raster(const grid_t &dims, const O* src_data):
            m_extent(dims),
            m_size(detail::array_product<M>(PointType<index_t,M>::iterator(dims)) * N),
            m_data(new O[m_size]),
            m_domain((I)0, (coord_t)dims - (coord_t)1) {
        std::copy(src_data, src_data + m_size, m_data.get());
    }
    
    /** 
     * Construct a Raster with `dims` elements along each axis.
     * Fill the raster with the data in `src_data`. `src_data` must have enough 
     * elements to fill the raster completely; in other words, the length of `src_data`
     * must be the product of the elements of `dims`. 
     * 
     * @param dims Number of samples along each axis.
     * @param src_data Data to copy into this raster.
     * @param domain Desired boundary of the data region. The lower and upper 
     * coordinates correspond to the exact coordinates of the most extreme data 
     * points along each axis.
     */
    Raster(const grid_t &dims, const O* src_data, const Rect<I,M> &domain):
            m_extent(dims),
            m_size(detail::array_product<M>(PointType<index_t,M>::iterator(dims)) * N),
            m_data(new O[m_size]),
            m_domain(domain) {
        std::copy(src_data, src_data + m_size, m_data.get());
    }
    
    ////////// Methods //////////
    
    /**
     * Set the sample value to be used for "out-of-bounds" samples when 
     * `EDGE_CONSTANT` behavior is used.
     */
    void setAbyss(sample_t val) {
        m_abyss = val;
    }
    
    /**
     * Set the value of the raster at the grid coordinate `idx`. 
     * @param idx Coordinate of datapoint to set.
     * @param val New value of datapoint.
     */
    void set(const grid_t &idx, const sample_t &val) {
        if (contains_gridpt(idx)) {
            O *p = m_data.get() + this->template index<EDGE_CONSTANT>(idx);
            //xxx no worky with dynamic
            *(sample_t*)p = val;
        }
    }
    
    /**
     * Copy the data in `region` into the array `dest`, in row-major 
     * (first dimension consecutive) order. Templated over edge sampling
     * behavior.
     * 
     * @tparam Edge  Edge sampling behavior.
     * @param dest   Destination array.
     * @param region Grid region of this Raster to copy from.
     */
    template <EdgeBehavior Edge>
    void copy(sample_t *dest, const Rect<index_t,M> &region) const {
        GridIterator<index_t,M> i = GridIterator<index_t,M>(region);
        GridIterator<index_t,M> end = i.end();
        for (; i != end; ++i, ++dest) {
            *dest = this->template sample_discrete<Edge>(*i);
        }
    }
    
    /**
     * Copy the data in `region` into the array `dest`, in row-major
     * (first dimension consecutive) order. Dynamic edge sampling behavior.
     * 
     * @param dest   Destination array.
     * @param region Grid region of this raster to copy from.
     * @param edge   Edge sampling behavior.
     */
    void copy(sample_t *dest, const Rect<index_t,M> &region, EdgeBehavior edge) const {
        GridIterator<index_t,M> i = GridIterator<index_t,M>(region);
        GridIterator<index_t,M> end = i.end();
        for (; i != end; ++i, ++dest) {
            *dest = this->sample_discrete(*i, edge);
        }
    }
    
    /**
     * Sample this Raster at `pt`.
     * 
     * @tparam Edge   Edge sampling behavior.
     * @tparam Interp Sample interpolation strategy.
     * @param pt      Point to sample.
     * @return Sampled data.
     */
    template <EdgeBehavior Edge, Interpolation Interp>
    inline sample_t sample(const coord_t &pt) const {
        return detail::_ImplRasterSample<I,O,M,N,Edge,Interp>::sample(this, toGridSpace(pt));
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
        pt = toGridSpace(pt);
        
        switch (interp) {
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
            return m_abyss;
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
        switch (Edge) {
            case EDGE_CONSTANT:
                if (contains_gridpt(pt)) {
                    offs = this->template index<EDGE_CONSTANT>(pt);
                } else {
                    return m_abyss;
                }
                break;
            default:
                offs = this->template index<Edge>(pt);
                break;
        }
        return PointType<O,N>::from_ptr(m_data.get() + offs);
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
        
        switch (edge) {
        case EDGE_CONSTANT:
            if (!contains_gridpt(pt)) return m_abyss;
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
            return m_abyss;
        }
        
        return sample_t(m_data.get() + offs);
    }
    
    /**
     * @return `true` if the grid location `idx` is within the bounds of this
     * raster and has data, `false` otherwise.
     */
    bool contains_gridpt(const grid_t &idx) const {
        typedef PointType<index_t,M> grid_info;
        
        const index_t *p_i = grid_info::iterator(idx);
        const index_t *e_i = grid_info::iterator(m_extent);
        for (index_t i = 0; i < M; i++) {
            index_t c = p_i[i];
            if (c < 0 or c >= e_i[i]) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Number of samples along each axis.
     */
    inline grid_t dataExtents() const {
        return m_extent;
    }
    
    /**
     * Boundary of the data region; upper and lower coordinates correspond exactly
     * to the coordinates of the most extreme sample positions along each axis.
     */
    inline Rect<I,M> domain() const {
        return m_domain;
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
    
    /**
     * Number of grid sample points in this raster.
     */
    inline index_t samplecount() const {
        return m_size;
    }
    
protected:
    
    inline coord_t toGridSpace(const coord_t &pt) const {
        return (coord_t)m_extent * m_domain.unmap(pt);
    }
    
    template <EdgeBehavior Edge>
    inline index_t index(const grid_t &c) const {
        // specialized where M=1
        return detail::_ImplRasterIndex<grid_t,M,Edge>::index(this->m_extent, c) * outputDimension();
    }
};

} // namespace geom

#endif /* RASTER_H_ */
