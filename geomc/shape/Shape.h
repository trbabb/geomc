#pragma once
/*
 * Shape.h
 *
 *  Created on: Oct 7, 2010
 *      Author: tbabb
 */

#include <geomc/shape/ShapeTypes.h>
#include <geomc/shape/Rect.h>

// todo: other concepts:
// - separate contains(pt) from sdf? ("Closed?")
// - could make RayIntersectable<...> into Intersectable<OtherShape>
//   - note that the return type for Ray is special, and might be for other shapes too
// todo: add a unified translation op that can work on most shapes
//   (natural for all but Frustum)

// todo: other compositional shapes?
// x Intersected<Shape1, Shape2>
//   - ray intersection is the intersection of ray intervals (RayIntersectable)
//   x convex_support is the point with lower dot with `d` (Convex)?
//      > nope, that's wrong
//      > if that were true, then you could trivially construct an intersection
//        shape and then check if its bbox is empty, instead of gjk
//      > imagine with two intersecting spheres— you get surface points
//        on a non-convex shape.
//   x sdf is max(sdf(a), sdf(b)) (SdfEvaluable)?
//      > nope, that's a bound, not a true sdf
//   x project() is the more inward of the two projection points (Projectable)?
//      > also wrong. consider two overlapping spheres; a lenticular shape.
//        many points should project to the ridge; that ridge point is usually
//        not the projection point on either of the two shapes.
//      > haven't proven this route impossible, though. something magic
//        might be possible by a combination of project() and support()
//        (but quite unlikely as it would be novel).
//      - you want: douglas-rachford algorithm
//        - https://x.com/gabrielpeyre/status/1791332821459472643
//   - some shapes are closed under intersection (Rect)
// - minkowski sum
//   - ray intersection is interval minkowski sum (implemented)
//   - convex_support() is the sum of the two pts
//   x I think project() is b.project(p - a.project(p)) + a.project(p)
//     x this is incorrect
//     > consider any shape which isn't centered on the projected point
//     > this is a very general inversion problem: find a point on the shape
//       at which to center the added shape, such that p's projection distance
//       is minimal. there is no general way to know, without searching, which of
//       all possible points in the shape that is.
//     > if this were true it would again suggest a trivial alternative to gjk
//       (check minkowski_shape.contains(origin))


namespace geom {

/** @addtogroup shape
 *  @{
 */


/**
 * @brief Base class describing shapes with finite extents in N dimensions.
 *
 * Uses the curiously-recurring template pattern to perform static polymorphism.
 * Override `bounds()` in the derived implementation.
 */
template <typename T, index_t N, typename Derived>
class Bounded : public Dimensional<T,N> {
    public:
    
    /**
     * @brief Produces an axis-aligned box completely enclosing this shape.
     */
    Rect<T,N> bounds() const {
        Derived* self = (Derived*) this;
        return self->bounds();
    }
};


/**
 * @brief Base class describing convex shapes in N-dimensional space.
 *
 * Uses the curiously-recurring template pattern to perform static polymorphism.
 * Override `convex_support()` in the derived implementation.
 */
template <typename T, index_t N, typename Derived>
class Convex : public Bounded<T,N,Derived> {
    public:
    
    using typename Dimensional<T,N>::point_t;
    
    /**
     * @brief Geometric convex support function.
     *
     * Returns the point on the surface of this convex shape that is furthest 
     * along direction `d` (i.e., has the highest dot product with `d`). 
     * 
     * All shapes which implement this function automatically support geometrical
     * intersection tests with any other Convex object.
     * 
     * @param d Direction along which to find a support plane.
     * @return A point on the surface of this convex shape.
     */
    point_t convex_support(point_t d) const {
        static_assert(&Derived::convex_support != &Convex<T,N,Derived>::convex_support,
            "Derived class must provide a convex_support() implementation");
        Derived* self = (Derived*) this;
        return self->convex_support(d);
    }
    
    /**
     * @brief Convex shape overlap test.
     * 
     * @return True if and only if this convex shape overlaps `other`; false otherwise.
     */
    template <ConvexObject Shape>
    requires (Shape::N == N) and std::same_as<T, typename Shape::elem_t>
    bool intersects(const Shape& other) const {
        return geom::intersects(
            as_any_convex(*this),
            as_any_convex(other)
        );
    }
    
    /**
     * @brief Produces an axis-aligned box completely enclosing this shape.
     *
     * The default implementation calls `convex_support()` along each of the 
     * principal axes to find the extents.
     */
    Rect<T,N> bounds() const {
        Rect<T,N> b;
        T* lo = PointType<T,N>::iterator(b.lo);
        T* hi = PointType<T,N>::iterator(b.hi);
        for (index_t i = 0; i < N; ++i) {
            point_t axis, x;
            const Derived* self = (const Derived*) this;
            PointType<T,N>::iterator(axis)[i] =  1;
            x = self->convex_support(axis);
            hi[i] = PointType<T,N>::iterator(x)[i];
            PointType<T,N>::iterator(axis)[i] = -1;
            x = self->convex_support(axis);
            lo[i] = PointType<T,N>::iterator(x)[i];
        }
        return b;
    }
    
};


/**
 * @brief Base class describing N-dimensional shapes which can be intersection-tested
 * with a Ray.
 * 
 * Uses the curiously-recurring template pattern to perform static polymorphism.
 */
template <typename T, index_t N, typename Derived>
class RayIntersectable {
public:
    /**
     * @brief Ray/shape intersection.
     *
     * Return the possibly-empty range of ray parameters `s` such that the Ray
     * overlaps the shape at `r.origin + s * r.direction`.
     * 
     * For non-convex objects with multiple regions of overlap, this method should
     * return the interval containing the smallest positive value.
     */
    Rect<T,1> intersect(const Ray<T,N>& r) const {
        static_assert(&Derived::intersect != &RayIntersectable<T,N,Derived>::intersect,
            "Derived class must provide an intersect(Ray) implementation");
        Derived* self = (Derived*) this;
        return self->intersect(r);
    }
};


/**
 * @brief Base class describing N-dimensional shapes which implement a signed
 * distance function.
 * 
 * Uses the curiously-recurring template pattern to perform static polymorphism.
 * Override `sdf()` in the derived implemenation.
 */ 
template <typename T, index_t N, typename Derived>
class SdfEvaluable {
public:
    
    /**
     * @brief Signed distance function.
     *
     * Compute a signed distance to the nearest surface point on the shape.
     * Points on the exterior have positive value; points on the interior have negative
     * value; surface points have sdf value 0.
     */
    inline T sdf(Vec<T,N> p) const {
        static_assert(&Derived::sdf != &SdfEvaluable<T,N,Derived>::sdf,
            "Derived class must provide an sdf() implementation");
        Derived* self = (Derived*) this;
        return self->sdf(p);
    }
    
    /**
     * @brief Shape-point overlap test.
     * 
     * Return `true` if the point `p` is on the surface or interior of the shape,
     * `false` otherwise.
     */
    inline bool contains(Vec<T,N> p) const {
        Derived* self = (Derived*) this;
        return self->sdf(p) <= 0;
    }
};


/**
 * @brief Base class describing N-dimensional shapes which implement the ability
 * to project an arbitrary point to the nearest point on their surface.
 *
 * Uses the curiously-recurring template pattern to perform static polymorphism.
 * Override `project()` and at least one of `sdf()` or `contains()` in the derived
 * implementation.
 */
template <typename T, index_t N, typename Derived>
class Projectable: public SdfEvaluable<T,N,Derived> {
public:
    
    /**
     * @brief Nearest point on the surface of the shape.
     *
     * Compute the point on the boundary of the shape which is closest to `p`. 
     */
    Vec<T,N> project(Vec<T,N> p) const {
        static_assert(&Derived::project != &Projectable<T,N,Derived>::project,
            "Derived class must provide a project() implementation");
        Derived* self = (Derived*) this;
        return self->project(p);
    }
    
    T sdf(Vec<T,N> p) const {
        static_assert(
            &Derived::sdf      != &Projectable<T,N,Derived>::sdf or
            &Derived::contains != &SdfEvaluable<T,N,Derived>::contains,
            "Derived class must override at least one of `sdf()` or `contains()`");
        Derived* self = (Derived*) this;
        return (self->contains(p) ? -1 : 1) * self->project(p).dist(p);
    }
    
    /**
     * @brief Unit-length outward-facing direction.
     *
     * For any point, return the direction which points directly away from the
     * nearest point on the shape, away from its interior.
     *
     * This should be the same as the gradient of the sdf().
     */
    Vec<T,N> normal(Vec<T,N> p) const {
        Derived* self = (Derived*) this;
        Vec<T,N> n = (p - self->project(p)).unit();
        return (self->contains(p) ? -1 : 1) * n;
    }
    
    /**
     * @brief Nearest point on the interior of the shape.
     * 
     * If `p` is on the interior of the shape, return `p` unaltered; otherwise
     * orthogonally project `p` to the shape's surface.
     */
    Vec<T,N> clip(Vec<T,N> p) const {
        Derived* self = (Derived*) this;
        return self->contains(p) ? p : project(p);
    }
};

/**
 * @brief Helper class which virtualizes the static polymorphism of Convex shapes.
 *
 * This allows the various Shape classes to interoperate with pointers in specific
 * cases where it's needed. An example is in the implementation of 
 * `gjk_intersect(Shape1, Shape2)`, which should not instantiate a template for each
 * combination of shape classes.
 *
 * To convert a shape to an AnyConvex, call the convenience function
 * `as_any_convex(Shape)`.
 */
template <typename T, index_t N>
class AnyConvex: public Convex<T,N,AnyConvex<T,N>> {
public:
    using typename Convex<T,N,AnyConvex<T,N>>::point_t;
    
    point_t convex_support(point_t d) const {
        return _impl_convex_support(d);
    }
    
    Rect<T,N> bounds() const {
        return _impl_bounds();
    }
protected:
    virtual point_t   _impl_convex_support(Vec<T,N> p) const = 0;
    virtual Rect<T,N> _impl_bounds() const = 0;
};

/// Implementation of AnyConvex for a specific Shape.
template <ConvexObject Shape>
class AnyConvexImpl: public AnyConvex<typename Shape::elem_t, Shape::N> {
public:
    using typename Dimensional<typename Shape::elem_t, Shape::N>::elem_t;
    using Dimensional<typename Shape::elem_t, Shape::N>::N;
    using T = elem_t;
    using typename AnyConvex<T,N>::point_t;
    
    Shape shape;
    
    AnyConvexImpl(const Shape&  s):shape(s) {}
    AnyConvexImpl(const Shape&& s):shape(s) {}
    
protected:
    virtual point_t _impl_convex_support(Vec<T,N> p) const {
        return shape.convex_support(p);
    }
    
    virtual Rect<T,N> _impl_bounds() const {
        return shape.bounds();
    }
};

/**
 * @brief Wrap `s` in a virtual class which implements the Convex concept.
 */
template <typename Shape>
inline AnyConvexImpl<Shape> as_any_convex(const Shape& s) {
    return AnyConvexImpl<Shape>(s);
}

}
