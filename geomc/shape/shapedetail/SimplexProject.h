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
Vec<T,N> _project_to_face(
        const Simplex<T,N>& simplex,
        SimplexFace<T,N> face,
        Vec<T,N> point,
        Vec<T,N> to_p)
{
    Vec<T,N> out;
    // we pass this bc the facefinding algorithm already has it:
    // to_p = p - s.pts[best_face->included[best_n]];
    if (face.n == 0) {
        // projects to a vertex; the projection is the vertex itself.
        out = simplex[face.included[0]];
    } else if (face.n == 1) {
        // projects to an edge
        const Vec<T,N>& v = face.spanning()[0];
        out = to_p.project_on(v) + simplex.pts[face.included[1]];
    } else if (N > 2 and face.n == N - 1) {
        // p projects to a plane.
        // subtract the normal component of the vector from the face to the point
        out = point - to_p.project_on(face.normal());
    } else if (face.n == N) {
        // projects to the interior of the simplex; p unchanged
        out = point;
    } else if constexpr (N > 3) {
        // this case is only possible with N > 3.
        // project to the lowest-dimension subspace (nullspace or spanning basis).
        index_t n_null = face.n_null();
        bool proj_null = n_null < face.n;
        Vec<T,N>* proj_basis = proj_null ? face.nullspace() : face.spanning();
        index_t proj_n = proj_null ? n_null : face.n;
        if (simplex.n < N + 1 or not proj_null) {
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
        out += simplex.pts[face.included[face.n]];
    }
    
    return out;
}


template <typename T, index_t N>
struct ProjectionResult {
    // the candidate/resultant face to which the point projects
    SimplexFace<T,N> face;
    // intermediate candidate has been found?
    bool     exists = false;
    // the point on the face to which the input point projects
    Vec<T,N> projected_point;
    // squared distance from the point to its projection on the face
    T        distance2 = std::numeric_limits<T>::max();
    // in the process of search, has a degenerate face been encountered?
    bool     is_degenerate = false;
    // whether the point is inside a full simplex
    bool     contains = false;
    // whether the point is on the back side of a boundary face.
    // this can happen if we are not yet sure whether the point is inside
    // a full simplex
    bool     backfacing = false;
    
    operator bool() const {
        return exists;
    }
    
    void update(T d2, Vec<T,N> point, const SimplexFace<T,N>& face, bool new_backfacing) {
        if (not exists // no previous candidates
            // this candidate is closer
            or (d2 < distance2 and new_backfacing <= backfacing)
            // frontface always beats backface
            or (not new_backfacing and backfacing))
        {
            this->exists          = true;
            this->face            = face;
            this->projected_point = point;
            this->distance2       = d2;
            this->backfacing      = new_backfacing;
        }
    }
    
};


template <typename T, index_t N>
struct SimplexProjection {
    const Simplex<T,N>&   simplex;
    Vec<T,N>              p;
    ProjectionOp          op;
    SimplexFaces          check_faces;
    ProjectionResult<T,N> result;
    
    SimplexProjection(
            const Simplex<T,N>& s,
            Vec<T,N> p,
            ProjectionOp op,
            SimplexFaces check_faces=SimplexFaces::ALL): 
                simplex(s),
                p(p),
                op(op),
                check_faces(check_faces),
                result {
                    .projected_point = p
                }
    {
        _find_nearest_face(SimplexFace<T,N>{s, s.n > 2}, p - s.pts[s.n - 1], 0);
    }
    
    Vec<T,N> projected_point() const {
        return result.projected_point;
    }
    
    T distance() const {
        return std::sqrt(result.distance2);
    }
    
protected:
    
    // return whether the point is interior to the subface
    bool _find_nearest_face(const SimplexFace<T,N>& face, Vec<T,N> to_p, int depth) {
        bool all_inside = true;
        bool backfacing = depth > 0 and face.normal().dot(to_p) < 0;
        if (backfacing and op == ProjectionOp::PROJECT and depth == 1 and face.n == N - 1) {
            // we're on the interior side of a top-level boundary face, and we might project
            // the point to one such face if we're inside a full simplex. we don't need to check 
            // its subfaces, though, because you can't project to a subface (e.g., an edge)
            // from the interior of a simplex. we also shouldn't early exit because we're on
            // the backside! (in the case where we're clipping, interior faces are not
            // projection candidates).
            
            // because we don't boundary check in this case, we could in principle project
            // to the face's plane well outside the simplex. if we are truly inside
            // the simplex, then some later face will be closer and occlude that hit. but we
            // have to be sure that any frontfacing hit always clobbers a backfacing best-hit,
            // since the backfacing plane might skim closer to the point than the frontfacing
            // face it should project to.
        } else {
            if (backfacing) {
                // point is on the interior side of a boundary, 
                // so it does not project to this subface.
                return true;
            }
            if (face.n == 0) {
                // simplex is a point and has no subfaces to check
                result.update(to_p.mag2(), simplex[face.included[0]], face, false);
                return false;
            }
            
            SimplexFace<T,N> nearest_face;
            index_t n_to_check = face.n + (check_faces == SimplexFaces::SKIP_BACKFACE ? 0 : 1);
            for (index_t i = 0; i < n_to_check; ++i) {
                // check each subface for ownership of the projected point.
                SimplexFace<T,N> subface = face.exclude(i);
                
                // re-root p if we're checking the backface
                Vec<T,N> to_sub_p = (i == face.n) ? (to_p - face.bases[i - 1]) : to_p;
                if (subface.normal().mag2() == 0) [[unlikely]] {
                    result.is_degenerate = true;
                    continue; // degenerate subface
                }
                
                bool is_interior = _find_nearest_face(subface, to_sub_p, depth + 1);
            
                if (not is_interior) all_inside = false;
            }
        }
        
        if (all_inside) {
            // the point is on the interior of all our subfaces, therefore it projects to us.
            if (face.n < N) {
                Vec<T,N> projected_pt = _project_to_face(simplex, face, p, to_p);
                T d2 = p.dist2(projected_pt);
                result.update(d2, projected_pt, face, backfacing);
            } else {
                // if we're a full simplex, we don't project at all:
                //   - if projecting, we project to a top-level boundary face, not the root
                //   - if clipping, no change to the point at all.
                // we do know we're inside the simplex, though:
                result.contains = true;
                if (op == ProjectionOp::CLIP) {
                    result.update(0, p, face, true);
                }
            }
        }
        // this can be backfacing when we're checking the top-level boundary faces;
        // i.e. the point is interior to a 2D tri or a 3D tet
        return backfacing;
    }
    
public:
    
    /**
     * @brief Compute a direction normal to the projected face.
     * 
     * In some cases, there may be more numerically stable methods for obtaining
     * the normal than subtracting the projected point from the original point.
     */
    Vec<T,N> normal_direction() {
        Vec<T,N> out_normal = p - result.projected_point;
        index_t best_n = result.face.n; // dimension of the projection target simplex
        Vec<T,N> to_p = p - simplex.pts[result.face.included[best_n]];
        if (best_n == 1) {
            // projects to an edge
            const Vec<T,N>& v = result.face.spanning()[0];
            if constexpr (N == 3) {
                // we can use the cross product to get a direction to the edge.
                // numerically this is better than subtracting from the projected point
                // (catastrophic cancellation)
                Vec<T,N> c = v ^ to_p;
                out_normal = c ^ v;
            } else if constexpr (N == 2) {
                // there are only two orthogonal directions to a line in 2D
                Vec<T,N> d = v.right_perpendicular();
                out_normal = d.dot(to_p) >= 0 ? d : -d;
            }
        } else if (N > 2 and best_n == N - 1) {
            // p projects to a hyperplane. return its normal, which we explicitly know
            Vec<T,N> normal = result.face.nullspace()[0];
            out_normal = normal.dot(to_p) >= 0 ? normal : -normal;
        }
        
        return out_normal;
    }
    
    /**
     * @brief Extract the face to which `p` projects.
     */
    Simplex<T,N> projected_face() const {
        Simplex<T,N> out;
        for (index_t i = 0; i <= result.face.n; ++i) {
            out.pts[i] = simplex.pts[result.face.included[i]];
        }
        out.n = result.face.n + 1;
        return out;
    }

};

} // namespace detail
} // namespace geom
