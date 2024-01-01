#pragma once

#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Orthogonal.h>
#include <geomc/shape/Simplex.h>
#include <geomc/Storage.h>

#if DEBUG_INTERSECTION
#include <geomc/shape/shapedetail/GJKDebug.h>
#endif

namespace geom {

template <typename T, index_t N>
struct SimplexFace {
    
    // basis layout:
    // [spanning] [nullspace] [normal]
    // as the simplex shrinks, the nullspace grows. This layout ensures:
    // (1) that we don't have to copy the entire nullspace when we exclude a basis, and
    // (2) that the existing nullspace and the spanning basis are contiguous,
    // which means we can compute the nullspace by passing a single pointer.
    
    index_t  n = 0; // number of bases in the face (one less than the number of verts)
    Vec<T,N> bases[N];
    index_t  included[N + 1];
    
    SimplexFace() = default;
    
    SimplexFace(const SimplexFace<T,N>& other) = default;
    
    SimplexFace(const Simplex<T,N>& s, bool compute_nullspace=true):n(s.n - 1) {
        Vec<T,N>* basis = spanning();
        // the origin of the simplex is its final vertex:
        for (index_t i = 0; i < n; ++i) {
            basis[i] = s.pts[i] - s.pts[n];
            included[i] = i;
        }
        included[n] = n;
        if (n < N and compute_nullspace) {
            // compute the orthogonal complement
            geom::nullspace(basis, n, nullspace());
        }
    }
    
    index_t n_null() const {
        return N - n;
    }
    
    Vec<T,N>& normal() {
        return bases[N - 1];
    }
    
    const Vec<T,N>& normal() const {
        return bases[N - 1];
    }
    
    Vec<T,N>* spanning() {
        return bases;
    }
    
    const Vec<T,N>* spanning() const {
        return bases;
    }
    
    Vec<T,N>* nullspace() {
        return bases + n;
    }
    
    const Vec<T,N>* nullspace() const {
        return bases + n;
    }
    
    SimplexFace<T,N> exclude(index_t i) const {
        SimplexFace<T,N> out;
        out.n = n - 1;
        std::copy(bases,    bases    + N,     out.bases);
        std::copy(included, included + N + 1, out.included);
        // perform this shuffle:
        //  ejected <- final basis <- normal
        out.bases[i]     = bases[n - 1];
        out.bases[n - 1] = normal();
        // update the index map. note that the root point is always included
        // and must stay at the end, so we move the final non-root point into
        // the vacated hole, and then shift the root point down by one:
        out.included[i]     = included[n - 1]; // xxx handle exclude root
        out.included[n - 1] = included[n];
        if (i >= n) {
            // the excluded vertex is the origin of the basis.
            // re-root the basis at the next vertex.
            for (index_t k = 0; k < out.n; ++k) {
                out.bases[k] -= bases[included[out.n]];
            }
        }
        Vec<T,N> ejected = bases[i];
        if (out.n > 0) {
            // compute the normal to the new face. it's orthogonal to the new basis,
            // and also lies in the space of the parent face (orthogonal to the parent's
            // nullspace).
            out.normal() = orthogonal(out.bases);
            // make sure the normal faces outward:
            if (out.normal().dot(ejected) >= 0) {
                out.normal() = -out.normal();
            }
        } else {
            // the simplex is a point. we have only one choice left for the
            // new normal, and it's in the space of the parent simplex; i.e.
            // the basis we just ejected.
            out.normal() = -ejected; // xxx: make sure this is right in the case where we re-root
        }
        return out;
    }
    
};

/**
 * @brief Find the sub-face to which the point `p` projects. If a new best result
 * is found, assign it to `best_face`.
 * 
 * @param face The sub-face to check for projection.
 * @param best_face The current best face to which `p` projects.
 * @return `true` iff the point falls on the interior side of the face.
 */
template <typename T, index_t N>
bool find_nearest_face(
        SimplexFace<T,N> face,
        std::optional<SimplexFace<T,N>>* best_face,
        Vec<T,N> to_p,
        bool is_root,
        bool* is_degenerate=nullptr)
{
    // xxx: todo: handle excluding the final vertex
    bool frontfacing = is_root or face.normal().dot(to_p) >= 0;
    if (not frontfacing) return true;  // point is on the interior
    // early exit: a neighbor face already owns:
    if (best_face->has_value() and face.n <= (*best_face)->n) return false;
    size_t n_containing_faces = 0;
    bool all_inside = true;
    for (index_t i = 0; i < face.n; ++i) {
        // check each subface for ownership of the projected point.
        SimplexFace<T,N> subface = face.exclude(i);
        if (subface.normal().mag2() == 0) {
            if (is_degenerate) *is_degenerate = true;
            continue; // degenerate subface
        }
        if (find_nearest_face(subface, best_face, to_p, false)) {
            n_containing_faces += 1;
        } else {
            all_inside = false;
        }
    }
    if (all_inside and (not best_face->has_value() or n_containing_faces > (*best_face)->n)) {
        // the point is on the interior of all our subfaces, but exterior to
        // the face itself. therefore it projects to us.
        // we also check again against the dimension of the best face, in case a degenerate
        // subface reduced our dimension
        *best_face = face;
    }
    return false;
}


template <typename T, index_t N>
Vec<T,N> project_to_simplex(
        const Simplex<T,N>& s,
        Vec<T,N> p,
        Simplex<T,N>* out_simplex,
        bool* is_degenerate=nullptr)
{
    if (is_degenerate) *is_degenerate = false;
    std::optional<SimplexFace<T,N>> best_face = std::nullopt;
    Vec<T,N> to_p = p - s.pts[s.n - 1];
    find_nearest_face(SimplexFace<T,N>{s, s.n > 2}, &best_face, to_p, true, is_degenerate);
    
    Vec<T,N> out;
    index_t best_n = best_face->n;
    if (best_n == 0) {
        // projects to a vertex; the projection is the vertex itself.
        out = s[best_face->included[0]];
    } else if (best_n == 1) {
        // projects to an edge
        const Vec<T,N>& v = best_face->spanning()[0];
        out = to_p.project_on(v) + s.pts[best_face->included[1]];
    } else if (N > 2 and best_n == N - 1) {
        // p projects to a plane.
        // subtract the normal component of the vector from the face to the point
        out = p - to_p.project_on(best_face->normal());
    } else if (best_n == N) {
        // projects to the interior of the simplex; p unchanged
        out = p;
    } else if constexpr (N > 3) {
        // this case is only possible with N > 3.
        // project to the lowest-dimension subspace (nullspace or spanning basis).
        index_t n_null = best_face->n_null();
        bool proj_null = n_null < best_n;
        Vec<T,N>* proj_basis = proj_null ? best_face->nullspace() : best_face->spanning();
        index_t proj_n = proj_null ? n_null : best_n;
        if (s.n < N + 1 or not proj_null) {
            // the projection basis is not orthogonal if:
            // (a) the root simplex is not full, and we generated the nullspace
            //     on inititalization
            // (b) we're projecting to the simplex basis, which is in general
            //     not orthogonal
            orthogonalize(proj_basis, proj_n);
        }
        // project P onto the chosen basis
        for (index_t i = 0; i < proj_n; ++i) {
            out += to_p.project_on(proj_basis[i]);
        }
        if (proj_null) {
            // if we projected to the nullspace,
            // then subtract that component from the original vector,
            // leaving a vector in the space of the simplex.
            out = to_p - out;
        }
        // add back the origin of the simplex
        out += s.pts[best_face->included[best_n]];
    }
    
    if (out_simplex) {
        // copy the points corresponding to the selected bases
        // into the output simplex:
        for (index_t i = 0; i <= best_n; ++i) {
            out_simplex->pts[i] = s.pts[best_face->included[i]];
        }
        out_simplex->n = best_n + 1;
    }
    return out;
}

/**
 * @brief Compute the outward facing normal.
 * 
 * - more numerically stable for pts near the surface.
 *   - always ok for N <= 3
 *   - precision degrades when points are (nearly) coincident with an edge
 * - does not yet handle interior case
 * - does not yet handle final face case
 */
template <typename T, index_t N>
Vec<T,N> simplex_normal_direction(
        const Simplex<T,N>& s,
        Vec<T,N> p,
        Simplex<T,N>* out_simplex,
        bool* is_degenerate=nullptr)
{
    if (is_degenerate) *is_degenerate = false;
    std::optional<SimplexFace<T,N>> best_face = std::nullopt;
    Vec<T,N> to_p = p - s.pts[s.n - 1];
    find_nearest_face(SimplexFace<T,N>{s, s.n > 2}, &best_face, to_p, true, is_degenerate);
    
    Vec<T,N> out;
    index_t best_n = best_face->n; // dimension of the projection target simplex
    if (best_n == 0) {
        // projects to a vertex; the normal is the direction from that vertex to the point
        out = to_p;
    } else if (best_n == 1) {
        // projects to an edge
        const Vec<T,N>& v = best_face->spanning()[0];
        if constexpr (N == 3) {
            // we can use the cross product to get a direction to the edge
            Vec<T,N> c = v ^ to_p;
            out = c ^ v;
        } else if constexpr (N == 2) {
            // there are only two orthogonal directions to a line in 2D
            Vec<T,N> d = v.right_perpendicular();
            out = d.dot(to_p) >= 0 ? d : -d;
        } else {
            out = to_p - to_p.project_on(v);
        }
    } else if (N > 2 and best_n == N - 1) {
        // p projects to a hyperplane. return the normal
        Vec<T,N> normal = best_face->nullspace()[0];
        out = normal.dot(to_p) >= 0 ? normal : -normal;
    } else if (best_n == N) {
        // projects to the interior of the simplex; no normal
        out = {};
    } else if constexpr (N > 3) {
        // this case is only possible with N > 3.
        // project to the lowest-dimension subspace (nullspace or spanning basis).
        index_t n_null = best_face->n_null();
        bool proj_null = n_null < best_n;
        Vec<T,N>* proj_basis = proj_null ? best_face->nullspace() : best_face->spanning();
        index_t proj_n = proj_null ? n_null : best_n;
        if (s.n < N + 1 or not proj_null) {
            // the projection basis is not orthogonal if:
            // (a) the root simplex is not full, and we generated the nullspace
            //     on inititalization
            // (b) we're projecting to the simplex basis, which is in general
            //     not orthogonal
            orthogonalize(proj_basis, proj_n);
        }
        // project P onto the chosen basis
        for (index_t i = 0; i < proj_n; ++i) {
            out += to_p.project_on(proj_basis[i]);
        }
        if (not proj_null) {
            // if we projected onto the simplex basis,
            // compute the null direction
            out = to_p - out;
        }
    }
    
    if (out_simplex) {
        // copy the points corresponding to the selected bases
        // into the output simplex:
        for (index_t i = 0; i <= best_n; ++i) {
            out_simplex->pts[i] = s.pts[best_face->included[i]];
        }
        out_simplex->n = best_n + 1;
    }
    return out;
}


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
            d = simplex_normal_direction(*cur_simplex, {}, next_simplex, &was_degenerate);
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
