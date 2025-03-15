#pragma once

/*
 * ShapeTypes.h
 *
 *  Created on: Nov 11, 2010
 *      Author: tbabb
 */

#include <geomc/linalg/LinalgTypes.h>

/** 
 * @defgroup shape Shape
 * @brief Shape-related functions and classes.
 */
 
 /**
  * @defgroup spline Spline
  * @brief Splines and curves.
  * @ingroup shape
  */

namespace geom {
    
/**
 * @addtogroup shape 
 * @{
 */

/** 
 * @brief Array traversal order, specified in terms of which axes to increment first.
 */
enum ArrayOrder {
    /**
        * @brief Array traversal by incrementing the first coordinate in the 
        * innermost loop. 
        * 
        * For matrices, whose coordinates are ordered `(row, col)`, this 
        * represens column-major (i.e. "fortran") order. If coordinates
        * are `(x, y, ...)`, then this is row-major order.
        */
    ARRAYORDER_FIRST_DIM_CONSECUTIVE,
    /**
        * @brief Array traversal by incrementing the last coordinate in the
        * innermost loop.
        * 
        * For matrices, whose coordinates are ordered `(row, col)`, this
        * represents row-major (i.e. "C") order. If coordinates are `(x, y, ...)`,
        * then this is column-major order.
        */
    ARRAYORDER_LAST_DIM_CONSECUTIVE
};

// base shape concepts
template <typename T, index_t N, typename Derived> class Bounded;
template <typename T, index_t N, typename Derived> class Convex;
template <typename T, index_t N, typename Derived> class Projectable;
template <typename T, index_t N, typename Derived> class SdfEvaluable;
template <typename T, index_t N, typename Derived> class RayIntersectable;

// shapes
template <typename T, index_t N> class Plane;
template <typename T, index_t N> class Sphere;
template <typename T, index_t N> class Rect;
template <typename T, index_t N> class Cylinder;
template <typename T, index_t N> class Simplex;
template <typename T, index_t N> class Capsule;
template <typename T, index_t N> class SphericalCap;
template <typename Shape>        class Extruded;
template <typename Shape>        class Transformed;
template <typename Shape>        class Similar;
template <typename Shape>        class Frustum;
template <typename Shape>        class Dilated;
template <typename Shape>        class Hollow;

/// @} // end group shape

/**
 * @addtogroup spline
 * @{
 */

// splines
template <typename T, index_t N, typename Derived> class CubicSpline;
template <typename T, index_t N> class BezierSpline;
template <typename T, index_t N> class BSpline;
template <typename T, index_t N> class CatromSpline;
template <typename T, index_t N> class HermiteSpline;
template <typename T, index_t N> class PolynomialSpline;

// paths
template <typename T, index_t N> class BSplinePath;
template <typename T, index_t N> class CatromSplinePath;
template <typename T, index_t N> class BezierPath;
template <typename T, index_t N> class HermitePath;

/// @} // end group spline

/**
 * @addtogroup shape
 * @{
 */

/**
 * @brief Convenience typedef for transformed Rects.
 * @related Transformed
 */
template <typename T, index_t N>
using AffineBox = Transformed<Rect<T,N>>;

/**
 * @brief Convenience typedef for arbitrarily-oriented Rects.
 * @related Similar
 */
template <typename T, index_t N>
using Box = Similar<Rect<T,N>>;

/**
 * @brief Convenience typedef for a spherical shell.
 */
template <typename T, index_t N>
using SphericalShell = Dilated<Hollow<Sphere<T,N>>>;

/**
 * @brief A 2D circle.
 * 
 * @tparam T Coordinate type
 * 
 * @related Sphere
 */
template <typename T>
using Circle = Sphere<T,2>;

template <typename T, index_t N, typename Object, typename NodeData=void*>
    class KDTree;
template <typename T, index_t N, ArrayOrder Order=ARRAYORDER_FIRST_DIM_CONSECUTIVE>
    class GridIterator;


///////////// Shape Concepts //////////////

/**
 * @brief A shape which can efficiently report its bounding box.
 */
template <typename Shape>
concept BoundedObject = DimensionalObject<Shape> and requires (const Shape s) {
    { s.bounds() } -> std::convertible_to<Rect<typename Shape::elem_t, Shape::N>>;
};

/**
 * @brief A convex shape which implements its convex support function.
 *
 * ConvexObject implies BoundedObject, because the support function can be used
 * to compute the bounding box.
 */
template <typename Shape>
concept ConvexObject = DimensionalObject<Shape> and BoundedObject<Shape> and 
    requires (const Shape s)
{
    { s.convex_support(Vec<typename Shape::elem_t, Shape::N>()) } -> std::convertible_to<typename Shape::point_t>;
};

/**
 * @brief A shape which can be intersected by a ray.
 */
template <typename Shape>
concept RayIntersectableObject = DimensionalObject<Shape> and requires (const Shape s) {
    { s.intersect(Ray<typename Shape::elem_t, Shape::N>()) } -> std::convertible_to<Rect<typename Shape::elem_t, 1>>;
};

/**
 * @brief A shape which can test for point containment.
 */
template <typename Shape>
concept RegionObject = DimensionalObject<Shape> and requires (const Shape s) {
    { s.contains(typename Shape::point_t()) } -> std::convertible_to<bool>;
};

/**
 * @brief A shape which can report its exact signed distance field.
 *
 * SdfObject implies RegionObject, because the signed distance field can be used
 * to test for containment.
 */
template <typename Shape>
concept SdfObject = DimensionalObject<Shape> and RegionObject<Shape> and
    requires (Shape s, typename Shape::point_t p)
{
    { s.sdf(p) } -> std::convertible_to<typename Shape::elem_t>;
};

/**
 * @brief A shape which can measure its interior region.
 * 
 * For a 2D shape, this is the area; for a 3D shape, this is the volume;
 * higher-dimensional shapes measure their hypervolumes.
 *
 * These methods avoid the use of "area" and "volume" in their names to avoid
 * confusion with the `measure_boundary()` method. ("Area" would measure the interior
 * of a 2D shape, but the boundary of a 3D shape).
 */
template <typename Shape>
concept InteriorMeasurableObject = DimensionalObject<Shape> and
    requires (const Shape& s)
{
    { s.measure_interior() } -> std::convertible_to<typename Shape::elem_t>;
};

/**
 * @brief A shape which can measure its boundary.
 *
 * For a 2D shape, this is its perimeter; for a 3D shape, this is its surface area;
 * higher-dimensional shapes measure their (N-1)-hypervolumes.
 */
template <typename Shape>
concept BoundaryMeasurableObject = DimensionalObject<Shape> and
    requires (const Shape& s)
{
    { s.measure_boundary() } -> std::convertible_to<typename Shape::elem_t>;
};

/**
 * @brief A shape which can measure both its surface area and volume.
 */
template <typename Shape>
concept MeasurableObject = InteriorMeasurableObject<Shape> and BoundaryMeasurableObject<Shape>;

/**
 * @brief A shape which can perform an orthogonal projection to its boundary.
 *
 * `ProjectionObject`s can also report their normal vectors— the direction
 * away from the surface at any given point— and can implelement a clip(),
 * which projects a point onto the surface if and only if it is outside.
 *
 * ProjectionObject implies SdfObject, because the projection function
 * can be used to implement the signed distance field.
 */
template <typename Shape>
concept ProjectableObject = DimensionalObject<Shape> and SdfObject<Shape> and
    requires (const Shape s, typename Shape::point_t p)
{
    { s.project(p) } -> std::convertible_to<typename Shape::point_t>;
    { s.normal(p)  } -> std::convertible_to<typename Shape::point_t>;
    { s.clip(p)    } -> std::convertible_to<typename Shape::point_t>;
};


///////////// Intersection helpers //////////////

// these make reference to the Concepts and so need to be down here

// virtualization helper
template <typename T, index_t N> class AnyConvex;
template <ConvexObject Shape>    class AnyConvexImpl;

template <typename T, index_t N>
bool intersects(
    const AnyConvex<T,N>& shape_a,
    const AnyConvex<T,N>& shape_b
);

/// @} // addtogroup shape
    
} // namespace geom
