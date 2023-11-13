#pragma once

/*
 * geomc_defs.h
 *
 *  Created on: Dec 24, 2010
 *      Author: tbabb
 */

#include <cstddef>
#include <type_traits>
#include <limits>
#include <algorithm>
#include <stdint.h>

/**
 * @mainpage GEOMC Library
 * 
 * Geomc is a C++ template library for geometry and basic linear algebra. It is 
 * built to provide building blocks for 2D and 3D applications, though generalization
 * to `N` dimensions is a major design philosophy. 
 * 
 * Wherever possible, types are designed to interoperate intuitively via C++
 * arithmetic operators. Performance, templatization over element types and 
 * dimension, and minimization of dynamic memory allocations are all emphasized.
 * 
 * Source code
 * ===========
 * 
 * Download the geomc headers and source from:
 * http://github.com/trbabb/geomc
 * 
 */

// Shall a matrix check whether index is out-of-bounds on every access?

  /* #define GEOMC_MTX_CHECK_BOUNDS */

/**
 * Shall matrices verify correct dimensions before performing inter-matrix
 * operations? (recommended)
 */

#define GEOMC_MTX_CHECK_DIMS

/**
 * Shall matrices check and handle storage aliasing among matrices in
 * inter-matrix operations?
 */

#ifndef GEOMC_MTX_CHECK_ALIASING
#define GEOMC_MTX_CHECK_ALIASING
#endif

// Shall vectors include functions for outputting to streams?

#ifndef GEOMC_USE_STREAMS
#define GEOMC_USE_STREAMS 1
#endif

// Shall the <function> module include functions for outputting to streams?

#define GEOMC_FUNCTION_USE_STREAMS


#define PI  (3.141592653589793238462643383)
#define TAU (6.283185307179586476925286767)

#define DYNAMIC_DIM (0)

#define M_CLAMP(v,lo,hi) std::min(std::max((v),(lo)),(hi))

#define M_ENABLE_IF(cond) \
    typename std::enable_if<(cond), int>::type DUMMY=0

#define DERIVED_TYPE(base,derived) \
    typename std::enable_if< std::is_base_of< (base), (derived) >, (derived)>::type

#define REQUIRE_INHERIT(base,derived) \
    typename std::enable_if< std::is_base_of< (base), (derived) >, int>::type dummy=0

#define REQUIRE_INHERIT_T(base,derived) \
    typename std::enable_if< std::is_base_of<  base,   derived  >, int>::type

typedef std::ptrdiff_t index_t;


/** @brief Namespace of all `geomc` functions and classes. */
namespace geom {
    
    // storage fwd decls
    template <typename T, index_t N> struct Storage;
    template <typename T, index_t N> struct SizedStorage;
    template <typename T, index_t N> struct UnmanagedStorage;

};

#ifdef PARSING_DOXYGEN

/** @brief Functions to extend support of stdlib to geomc classes. */
namespace std { };

#endif
