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
 * @defgroup function
 * @brief Tools for constructing and sampling continuous-valued objects.
 */

// classes
template <typename I, typename O, index_t M, index_t N> class Raster;
template <typename I, typename O> class Expr;
template <typename T, index_t N>  class PerlinNoise;
template <typename T, index_t N>  class Path;

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

// perlin
typedef PerlinNoise<double,2> PerlinNoise2d;
typedef PerlinNoise<double,3> PerlinNoise3d;
typedef PerlinNoise<double,4> PerlinNoise4d;

typedef PerlinNoise<float,2> PerlinNoise2f;
typedef PerlinNoise<float,3> PerlinNoise3f;
typedef PerlinNoise<float,4> PerlinNoise4f;

// paths
typedef Path<double,2> Path2d;
typedef Path<double,3> Path3d;
typedef Path<double,4> Path4d;

typedef Path<float,2> Path2f;
typedef Path<float,2> Path3f;
typedef Path<float,2> Path4f;

// 1D functions
typedef Raster<float,float,1,1>   Curve1f;
typedef Raster<double,double,1,1> Curve1d;

// 2D images (float, float) -> (float, ...)
typedef Raster<float,float,2,1>   Imagef1f;
typedef Raster<float,float,2,2>   Imagef2f;
typedef Raster<float,float,2,3>   Imagef3f;
typedef Raster<float,float,2,4>   Imagef4f;

typedef Raster<double,double,2,1> Imaged1d;
typedef Raster<double,double,2,2> Imaged2d;
typedef Raster<double,double,2,3> Imaged3d;
typedef Raster<double,double,2,4> Imaged4d;

// 3D volumes (float, float, float) -> (float, ...)
typedef Raster<float,float,3,1>   Volumef1f;
typedef Raster<float,float,3,2>   Volumef2f;
typedef Raster<float,float,3,3>   Volumef3f;
typedef Raster<float,float,3,4>   Volumef4f;

typedef Raster<double,double,3,1> Volumed1d;
typedef Raster<double,double,3,2> Volumed2d;
typedef Raster<double,double,3,3> Volumed3d;
typedef Raster<double,double,3,4> Volumed4d;

/// @} // addtogroup function

} // namespace geom

#endif /* FUNCTIONTYPES_H_ */
