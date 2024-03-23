/*
 * FunctionTypes.h
 *
 *  Created on: Mar 16, 2013
 *      Author: tbabb
 */

#ifndef FUNCTIONTYPES_H_
#define FUNCTIONTYPES_H_

#include <geomc/linalg/LinalgTypes.h>

namespace geom {
    
/**
 * @defgroup function Function
 * @brief Tools for constructing and sampling continuous-valued objects.
 */

// classes
template <typename I, typename O, index_t M, index_t N> class Raster;
template <typename T, index_t N>     class PerlinNoise;
template <typename T, index_t N>     class Path;
template <typename T, index_t Bands> class SphericalHarmonics;
template <typename T, index_t Bands> class ZonalHarmonics;

/**
 * @addtogroup function
 * @{
 */

/// Raster edge-sampling behavior.
enum EdgeBehavior {
    /// Clip sample coordinates to the sampleable area, thus repeating edge values beyond the boundary
    EDGE_CLAMP,
    /// Wrap around sample coordinates to the opposite edge, thus tiling the sampled data beyond the boundary
    EDGE_PERIODIC,
    /// Mirror the sample coordinates across edges
    EDGE_MIRROR,
    /// Regions outside the sampleable area have a uniform, defined value (zero by default).
    EDGE_CONSTANT
};

/// Behavior for sampling between data points.
enum Interpolation {
    /// Return the data nearest to the sample point.
    INTERP_NEAREST,
    /// Linearly interpolate the nearest 2<sup>n</sup> data points.
    INTERP_LINEAR,
    /// Cubically interpolate the nearest 4<sup>n</sup> data points. 
    INTERP_CUBIC
};


/// @} // addtogroup function

} // namespace geom

#endif /* FUNCTIONTYPES_H_ */
