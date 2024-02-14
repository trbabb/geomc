#pragma once

#include <geomc/linalg/Vec.h>
#include <geomc/shape/ShapeTypes.h>

namespace geom {
namespace detail {

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
        std::copy(included, included + N + 1, out.included);
        Vec<T,N> ejected;
        if (i < n) {
            std::copy(bases, bases + N, out.bases);
            // perform this shuffle:
            //  ejected <- final basis <- normal
            out.bases[i]     = bases[n - 1];
            out.bases[n - 1] = normal();
            // update the index map. note that the root point is always included
            // and must stay at the end, so we move the final non-root point into
            // the vacated hole, and then shift the root point down by one:
            out.included[i]     = included[n - 1];
            out.included[n - 1] = included[n];
            ejected = bases[i];
        } else {
            // the excluded vertex is the current origin of the basis.
            // the new origin is the previous final basis; re-root all the axes
            for (index_t k = 0; k < out.n; ++k) {
                out.bases[k] = bases[k] - bases[n - 1];
            }
            // copy the nullspace
            std::copy(bases + n, bases + N, out.bases + out.n);
            // the ejected basis points from the new origin toward the old origin
            ejected = -bases[n - 1];
        }
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
            out.normal() = -ejected;
        }
        return out;
    }
    
};


enum struct ProjectionOp {
    CLIP,
    PROJECT
};

enum struct SimplexFaces {
    ALL,
    SKIP_BACKFACE
};


template <typename T, index_t N>
bool find_nearest_face(
        SimplexFace<T,N> face,
        std::optional<SimplexFace<T,N>>* best_face,
        Vec<T,N> to_p,
        ProjectionOp op,
        SimplexFaces check_faces,
        bool is_root,
        bool* is_degenerate=nullptr)
{
    bool frontfacing = is_root or face.normal().dot(to_p) >= 0;
    if (not frontfacing) return true;  // point is on the interior
    // early exit: a neighbor face already owns:
    if (best_face->has_value() and face.n < (*best_face)->n) return false;
    if (face.n == 0) {
        // simplex is a point and has no subfaces to check
        if (frontfacing) *best_face = face;
        return not frontfacing;
    }
    
    SimplexFace<T,N> nearest_face;
    T min_d2 = std::numeric_limits<T>::max();
    bool all_inside = true;
    index_t n_to_check = face.n + (check_faces == SimplexFaces::SKIP_BACKFACE ? 0 : 1);
    for (index_t i = 0; i < n_to_check; ++i) {
        // check each subface for ownership of the projected point.
        SimplexFace<T,N> subface = face.exclude(i);
        
        // re-root p if we're checking the backface
        if (i == face.n) to_p = to_p - face.bases[i - 1];
        if (subface.normal().mag2() == 0) {
            if (is_degenerate) *is_degenerate = true;
            continue; // degenerate subface
        }
        
        bool is_interior = find_nearest_face(
            subface,
            best_face,
            to_p,
            op,
            check_faces,
            false,
            is_degenerate
        );
        
        if (is_interior) {
            if (op == ProjectionOp::PROJECT and is_root and all_inside) {
                // keep track of the nearest wall, as long as we might be fully contained
                T d2 = to_p.project_on(subface.normal()).mag2();
                if (d2 < min_d2) {
                    nearest_face = subface;
                    min_d2 = d2;
                }
            }
        } else {
            all_inside = false;
        }
    }
    if (all_inside) {
        // the point is on the interior of all our subfaces, but exterior to
        // the face itself. therefore it projects to us.
        // we also check again against the dimension of the best face, in case a degenerate
        // subface reduced our dimension
        if ((is_root or min_d2 == 0) and op == ProjectionOp::PROJECT) {
            // we're fully inside the root simplex, and we care about which
            // wall we're closest to
            *best_face = nearest_face;
        } else {
            *best_face = face;
        }
    }
    return false;
}


template <typename T, index_t N>
Vec<T,N> project_to_simplex(
        const Simplex<T,N>& s,
        Vec<T,N> p,
        Simplex<T,N>* out_simplex,
        ProjectionOp op,
        bool* is_degenerate=nullptr)
{
    if (is_degenerate) *is_degenerate = false;
    std::optional<SimplexFace<T,N>> best_face = std::nullopt;
    Vec<T,N> to_p = p - s.pts[s.n - 1];
    
    find_nearest_face(
        SimplexFace<T,N>{s, s.n > 2}, 
        &best_face, 
        to_p,
        op,
        SimplexFaces::ALL,
        true, // is root
        is_degenerate
    );
    
    Vec<T,N> out;
    index_t best_n = best_face->n;
    to_p = p - s.pts[best_face->included[best_n]];
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
        ProjectionOp op,
        SimplexFaces check_faces,
        bool* is_degenerate=nullptr)
{
    if (is_degenerate) *is_degenerate = false;
    std::optional<SimplexFace<T,N>> best_face = std::nullopt;
    Vec<T,N> to_p = p - s.pts[s.n - 1];
    
    find_nearest_face(
        SimplexFace<T,N>{s, s.n > 2},
        &best_face,
        to_p,
        op,
        check_faces,
        true,
        is_degenerate
    );
    
    Vec<T,N> out;
    index_t best_n = best_face->n; // dimension of the projection target simplex
    to_p = p - s.pts[best_face->included[best_n]];
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


} // namespace detail
} // namespace geom
