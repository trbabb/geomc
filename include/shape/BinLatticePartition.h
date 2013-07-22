/*
 * BinLatticePartition.h
 *
 * Implements an N-dimensional, gridded, non-hierarchical spatial partition of cells.
 * Put another way, objects are kept in a grid of "buckets".
 *
 * Good for promximity/neighbor searches, when the density of objects
 * is roughly known.
 *
 * Spatial searches over many bins are likely to be inefficient, especially
 * with larger dimensions, as every bin within the search region is queried
 * regardless of whether it has contents.
 *
 * Bins are sparsely allocated and there is no memory cost for having
 * empty bins, and the structure itself does not have limited bounds.
 *
 * Multiple/duplicate objects may be placed at the same location.
 *
 * Plays nice in one dimension; point_t will revert to T instead of Vec<T,1>
 *
 *  Created on: Aug 19, 2012
 *      Author: tbabb
 */

// TODO: removing an object by value is difficult. (where did you stash it?)
//       maintain an inverse mapping?
//       - let's say that's up to the client, since
//         maintaining a mapping may not be necessary if, for example, the
//         location is directly calculable from the the stashed object.
// TODO: add current_bin() to ordinary iterator, and kill binsBegin()/binsEnd().
// TODO: implement move(item_t old_loc, point_t new_loc)? ... and move(iterator old_loc, point_t new_loc)?
// TODO: can we keep running tabs on the bounds of this object? (watch out for double-placed items)
//       no, probs not without a much smaerter (ordered) data structure.
//       foreach item you remove, where is the next boundary?
// TODO: add a filter iterator that selects by Area/Volume
// TODO: can we be more efficient for large region queries?
// TODO: can we have a place for a function which automatically maps <O> to point_t?
// TODO: equal_range? plus backward iterator?

////// other grid partition types //////

// TODO: have a region-based grid partition extensible to all Bounded.
//       if geo::intersects(obj, Rect<T,N>) is implemented, then this shall be called
//       for all bins that touch the box given by obj.bounds(). Otherwise, the object
//       shall be added to all bins touching obj.bounds(). The simpler case should be used
//       when obj is of type Rect.
// TODO: have a Vec-based grid partition that does not redundantly store location keys.
// TODO: take advantage of GridIterator
// TODO: move stuff to detail namespace.

// curious case: iType will not cover the same dynamic range that T will. what happens on overflow?
//               I think items can still be found; they will just fold into lower-order cells.

#ifndef BINLATTICEPARTITION_H_
#define BINLATTICEPARTITION_H_

#include <tr1/unordered_map>
#include <limits>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "linalg/Vec.h"
#include "shape/Rect.h"

namespace geom {

template<typename T, index_t N, typename O, typename iType = index_t>
class BinLatticePartition {

public:
    typedef typename Rect<iType,N>::point_t  bin_t;
    typedef typename Rect<T,N>::point_t      point_t;
    typedef std::pair< Vec<T,N>, O >         item_t;

protected:
    typedef std::tr1::unordered_multimap<bin_t, item_t> cellmap_t;

    // members
    cellmap_t _cells;
    T         _cellsize;

    // so we can transform the cellmap's iterator into something that returns pairs
    // of point:value and not bin:(point,value)
    static item_t& cell_unwrapper(typename cellmap_t::value_type &item){
        return item.second;
    }

    friend class region_iterator;
public:

    /*****************************
     * container iterators       *
     *****************************/

    // turn an iterator over bins into an iterator over location:value pairs
    // dark corners of c++ syntax here:
    //                        fnptr [  item_t  <-  pair(bin_t, item_t)  ]
    typedef boost::transform_iterator< item_t& (*)(typename cellmap_t::value_type&), typename cellmap_t::iterator> iterator;
    typedef typename cellmap_t::iterator bin_iterator;

    /*****************************
     * region iterator           *
     *****************************/

    class region_iterator: public boost::iterator_facade<region_iterator, item_t,
            boost::forward_traversal_tag> {
    public:
        /////// structors ///////

        // off-the-end iterator
        region_iterator() : grid(NULL),off_end(true) {}

        // iterator for a specific region
        explicit region_iterator(BinLatticePartition<T, N, O, iType> *grid,
                                 Rect<iType, N> region) :
                grid(grid),
                cellblock(region),
                off_end(false){
            cur_cell_iter = grid->_cells.find(region.min());
            if (cur_cell_iter == grid->_cells.end() and !gotoNextNonemptyCell(region.min())){
                off_end = true;
            }
        }

        operator iterator(){
            return iterator(cur_cell_iter, cell_unwrapper);
        }

        region_iterator end(){
            return region_iterator();
        }

        region_iterator begin(){
            if (grid == NULL){
                return region_iterator();
            } else {
                return region_iterator(grid, cellblock);
            }
        }
        
        bin_t current_bin(){
            if (not off_end){
                return cur_cell_iter->first;
            } else {
                return bin_t();
            }
        }
        
        Rect<iType,N> query_block(){
            return cellblock;
        }
        
    private:

        // TODO: understand invalidation (and its conventions).
        // TODO: const-ness?

        friend class boost::iterator_core_access;

        /////// methods ///////

        void increment() {
            if (!off_end){
                bin_t cur_cell = cur_cell_iter->first;
                cur_cell_iter++;
                // we must test whether we have hopped out of this bin; or the iterator will just keep
                // going on to another random one. which I find quite silly.
                if ((cur_cell_iter == grid->_cells.end() or cur_cell_iter->first != cur_cell) and !gotoNextNonemptyCell(cur_cell)){
                    off_end = true;
                }
            }
        }

        bool equal(region_iterator const & other) const {
            if (off_end and other.off_end){
                return true;
            } else {
                return other.cur_cell_iter == cur_cell_iter;
            }
        }

        item_t& dereference() const {
            return cur_cell_iter->second;
        }

        // iterate over cells in the region, checking each for existence.
        // return false if no nonempty cells exist within the region.
        bool gotoNextNonemptyCell(bin_t startCell) {
            bin_t thisCellCoords = startCell;
            while (true) {
                // increment cell within block
                bool cellWithinRegion = false;
                for (index_t inc_coord = 0; inc_coord < N; inc_coord++){
                    if (thisCellCoords[inc_coord] + 1 >= cellblock.max()[inc_coord]){
                        // don't exceeded maximum along this axis; can't increment.
                        // rollover this axis and attempt to increment the next.
                        thisCellCoords[inc_coord] = cellblock.min()[inc_coord];
                    } else {
                        thisCellCoords[inc_coord]++;
                        cellWithinRegion = true; // successfully incremented cell w/o going off edge
                        break;
                    }
                }

                // see if this cell has contents.
                if (cellWithinRegion){
                    typename cellmap_t::iterator i = grid->_cells.find(thisCellCoords);
                    if (i != grid->_cells.end()){
                        cur_cell_iter = i;
                        return true;
                    }
                } else {
                    break;
                }
            }
            return false;
        }

        /////// members ///////

        BinLatticePartition<T, N, O, iType> *grid; // iterator is off-end when this is null
        Rect<iType, N> cellblock; //query region
        typename cellmap_t::iterator cur_cell_iter; //iterator of currently inspected cell.
        bool off_end;

    }; //end nested class region_iterator

    /*****************************
     * structors                 *
     *****************************/

    explicit BinLatticePartition(T cellsize=1):_cellsize(cellsize) {};

    virtual ~BinLatticePartition() {};

    /*****************************
     * public methods            *
     *****************************/

    iterator begin(){
        return iterator(_cells.begin(), cell_unwrapper);
    }

    iterator end(){
        return iterator(_cells.end(), cell_unwrapper);
    }

    void insert(point_t location, O value) {
        bin_t bin = binOf(location);
        _cells.insert(std::pair<bin_t, item_t>(bin, item_t(location, value)));
    }

    void erase(const iterator &i){
        _cells.erase(i.base());
    }

    // delete an entire bin
    // TODO: WARNING: I do not think this does what you think it does.
    void erase(const bin_iterator &i){
        _cells.erase(i);
    }
    
    void clear(){
        _cells.clear();
    }

    region_iterator query(Rect<T,N> region){
        return region_iterator(*this, binRegion(region));
    }

    region_iterator query(point_t pt, T radius){
        point_t diag = point_t(radius);
        Rect<T,N> region = Rect<T,N>(pt-diag, pt+diag);
        Rect<iType,N> bins = binRegion(region);
        return region_iterator(this, bins);
    }

    size_t size(){
        return _cells.size();
    }

    //TODO: delete us
    bin_iterator binsBegin(){
        return _cells.begin();
    }

    bin_iterator binsEnd(){
        return _cells.end();
    }

    bin_t binOf(point_t location){
        return (bin_t)((location / _cellsize).floor());
    }

    protected:

    inline Rect<iType,N> binRegion(Rect<T,N> r){
        //high coordinate is, as always, one beyond the maximum cell to check
        return Rect<iType,N>(
                (bin_t)((r.min() / _cellsize).floor()),
                (bin_t)((r.max() / _cellsize).ceil())
        );
    }

};
//end class BinLatticePartition

};
//end namespace geom

#endif /* BINLATTICEPARTITION_H_ */
