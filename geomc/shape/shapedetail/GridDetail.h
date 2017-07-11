/*
 * GridDetail.h
 *
 *  Created on: Apr 8, 2013
 *      Author: tbabb
 */

#ifndef GRIDDETAIL_H_
#define GRIDDETAIL_H_

#include <boost/static_assert.hpp>
#include <geomc/shape/ShapeTypes.h>

namespace geom {
namespace detail {

template <ArrayOrder Order, index_t N>
struct _ImplArrayOrder {
    BOOST_STATIC_ASSERT(N != N);
};

// this redundant bullshit is necessary because c++ is badly designed.
// the partial specializations of the the static member variables each need to have
// their own corresponding class template, we can't just specialize the static
// variable by itself.
// at least we get the side benefit of explicitly disabling weird enum values for this class.

template <index_t N>
struct _ImplArrayOrder<ARRAYORDER_FIRST_DIM_CONSECUTIVE,N> {
    static const index_t dim_first;
    static const index_t dim_last;
    static const index_t dim_end;
    static const index_t dim_increment;
};

template <index_t N>
struct _ImplArrayOrder<ARRAYORDER_LAST_DIM_CONSECUTIVE,N> {
    static const index_t dim_first;
    static const index_t dim_last;
    static const index_t dim_end;
    static const index_t dim_increment;
};

// first dim "inner":

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_FIRST_DIM_CONSECUTIVE, N>::dim_first = 0;

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_FIRST_DIM_CONSECUTIVE, N>::dim_last = N-1;

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_FIRST_DIM_CONSECUTIVE, N>::dim_end = N;

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_FIRST_DIM_CONSECUTIVE, N>::dim_increment = 1;

// last dim "inner":

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_LAST_DIM_CONSECUTIVE, N>::dim_first = N-1;

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_LAST_DIM_CONSECUTIVE, N>::dim_last = 0;

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_LAST_DIM_CONSECUTIVE, N>::dim_end = -1;

template <index_t N>
const index_t _ImplArrayOrder<ARRAYORDER_LAST_DIM_CONSECUTIVE, N>::dim_increment = -1;

} // namespace detail
} // namespace geom

#endif /* GRIDDETAIL_H_ */
