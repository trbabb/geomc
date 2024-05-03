#pragma once

#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Orthogonal.h>
#include <geomc/shape/Simplex.h>
#include <geomc/Storage.h>

#if DEBUG_INTERSECTION
#include <geomc/shape/shapedetail/GJKDebug.h>
#endif

#include <geomc/shape/shapedetail/SimplexProject.h>


namespace geom {

template <typename T, index_t N>
struct Intersector {
    // afaict this is never reached in 2D or 3D
    // (2D, 3D, 4D, ... = 10, 20, 40, ...)
    index_t max_iterations = 10 * (1 << (N - 2));
    /// axis separating the shapes in the previous intersection test, if they overlapped.
    Vec<T,N> separation_axis  = Vec<T,N>::unit_x;
    /// number of GJK iterations taken to produce the last intersection result.
    index_t iterations        = 0;
    /// set to `true` iff the last intersection test resulted in a degenerate simplex
    bool was_degenerate       = false;

#if DEBUG_INTERSECTION
    FILE* debug_file = nullptr;
#endif

protected:
    // two "buffers":
    Simplex<T,N> simplex_a;
    Simplex<T,N> simplex_b;
    // front and back buffer:
    Simplex<T,N>* cur_simplex  = &simplex_a;
    Simplex<T,N>* next_simplex = &simplex_b;

public:

    const Simplex<T,N>& simplex() const { return *cur_simplex; }
    
    bool intersects(
        const AnyConvex<T,N>& shape_a,
        const AnyConvex<T,N>& shape_b)
    {
        // `a` is a point on the minkowski difference
        Vec<T,N> a = shape_a.convex_support( separation_axis) -
                     shape_b.convex_support(-separation_axis);
        // `d` is the previous search direction
        Vec<T,N> d = -a;
        // initialize the simplex with a single point:
        cur_simplex->n = 0;
        cur_simplex->insert(a);
        
        iterations = 0;
        while (true) {
            iterations += 1;
            a = shape_a.convex_support( d) - 
                shape_b.convex_support(-d);
            T k = a.dot(d);
            if (k < 0 or a == cur_simplex->pts[cur_simplex->n - 1]) {
                // we tried to search as far as we could in direction `d`,
                // but got no closer to the origin. 
                separation_axis = -a.project_on(d);
                // return whether the origin is inside the minkowski difference
                return k >= 0;
            }
            cur_simplex->insert(a);
            // project the origin onto the simplex,
            // putting the projection face into `next_simplex`
            detail::SimplexProjection<T,N> proj {
                *cur_simplex,
                {}, // origin
                detail::ProjectionOp::CLIP,
                detail::SimplexFaces::SKIP_BACKFACE
            };
            d = proj.normal_direction();
            *next_simplex  = proj.projected_face();
            was_degenerate = proj.result.is_degenerate;
            
#if DEBUG_INTERSECTION
            if (debug_file) {
                emit_splex_stage(debug_file, *cur_simplex, *next_simplex, a, -d);
            }
#endif
            std::swap(cur_simplex, next_simplex);
            if (cur_simplex->n == N + 1 or d.is_zero()) {
                // the simplex is full, and the origin is inside it.
                separation_axis = d;
                return true;
            }
            if (iterations > max_iterations or d.mag2() == 0) { return true; }
        }
    }

}; // struct Intersector


template <typename T, index_t N>
bool intersects(
    const AnyConvex<T,N>& shape_a,
    const AnyConvex<T,N>& shape_b)
{
    Intersector<T,N> intersector;
    return intersector.intersects(shape_a, shape_b);
}

} // namespace geom
