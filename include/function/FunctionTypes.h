/*
 * FunctionTypes.h
 *
 *  Created on: Mar 16, 2013
 *      Author: tbabb
 */

#ifndef FUNCTIONTYPES_H_
#define FUNCTIONTYPES_H_

#include "linalg/LinalgTypes.h"

namespace geom {

// classes
template <typename I, typename O, index_t N, index_t Channels> class Raster;
template <typename I, typename O> class Expr;
template <typename T, index_t N>  class PerlinNoise;
template <typename T, index_t N>  class Path;

enum EdgeBehavior {
    EDGE_CLAMP,
    EDGE_PERIODIC,
    EDGE_MIRROR,
    EDGE_CONSTANT
};

enum Interpolation {
    INTERP_NEAREST,
    INTERP_LINEAR,
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

}; // namespace geom

#endif /* FUNCTIONTYPES_H_ */
