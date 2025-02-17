#pragma once

#include <vector>
#include <geomc/shape/CubicSpline.h>

namespace geom {

// todo: allow closed vs open curves
// todo: standardize adding/removing knots,
//   so we can do things like wrap in arc-length parameterization,
//   distance-finding, ray intersection, etc.

/**
 * @addtogroup spline
 * @{
 */

/**
 * @brief Base class for a path defined by a sequence of concatenated splines.
 */
template <typename T, index_t N, CubicSplineObject<T,N> S, typename Derived>
class SplinePath : public Dimensional<T,N> {
    
          Derived& derived()       { return static_cast<      Derived&>(*this); }
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
    
public:
    
    using Spline = S;
    
    /**
     * @brief Return the index of the segment containing `s`.
     *
     * If there are no complete segments in the path, return `std::nullopt`.
     */
    std::optional<size_t> segment(T s) const {
        size_t n = derived().n_segments();
        if (n < 1) return {};
        return std::min(std::max<T>(s, 0), n - 1);
    }
    
    /**
     * @brief Evaluate the path at `s`.
     *
     * Segment `i` is evaluated for `s` in the range `[i, i + 1]`.
     * If `s` is outside the range `[0, n_segments()]`, the path is extrapolated
     * using the first or last segment.
     */
    VecType<T,N> operator()(T s) const {
        std::optional<size_t> i = segment(s);
        if (not i) { return {}; }
        Spline spline = derived()[*i].value();
        return spline(s - i);
    }
    
    /**
     * @brief Evaluate the derivative of the path at `s`, with respect to `s`.
     *
     * The range of `s` is the same as for `operator()`.
     */
    VecType<T,N> velocity(T s) const {
        std::optional<size_t> i = segment(s);
        if (not i) return {};
        Spline spline = (*this)[*i].value();
        return spline.velocity(s - i);
    }
    
    /**
     * @brief Evaluate the second derivative of the path at `s`, with respect to `s`.
     *
     * The range of `s` is the same as for `operator()`.
     */
    VecType<T,N> acceleration(T s) const {
        std::optional<size_t> i = segment(s);
        if (not i) return {};
        Spline spline = (*this)[*i].value();
        return spline.acceleration(s - i);
    }
    
    /**
     * @brief Compute the axis-aligned bounding box of the path.
     * 
     * O(n) on the number of segments.
     */
    Rect<T,N> bounds() const {
        Rect<T,N> r;
        for (size_t i = 0; i < derived().n_segments(); ++i) {
            r |= derived()[i].value().bounds();
        }
        return r;
    }
};


/**
 * @brief An extendable path defined by a sequence of guide knots.
 *
 * B-spline paths are concatenated cubic B-splines, where each segment shares
 * its first three control points with the previous segment.
 *
 * B-spline curves are C2 continuous, and do not necessarily pass through the
 * guide knots.
 */
template <typename T, index_t N>
class BSplinePath : public SplinePath<T,N,BSpline<T,N>,BSplinePath<T,N>> {
public:
    
    using Knot = VecType<T,N>;
    
    /// The sequence of guide knots.
    std::vector<VecType<T,N>> knots;
    
    /// Construct an empty path.
    BSplinePath() = default;
    /// Construct a path from a sequence of guide knots.
    BSplinePath(const std::initializer_list<VecType<T,N>> knots): 
        knots(knots) {}
    /// Construct a path having a single segment.
    BSplinePath(const BSpline<T,N>& segment):
        knots {segment.control_points(), segment.control_points() + 4}
    {}
    
    /**
     * @brief Return the number of complete segments in the path.
     *
     * A single segment is defined by four consecutive control points.
     */
    size_t n_segments() const {
        size_t n = knots.size();
        if (n < 4) return 0;
        return knots.size() - 3;
    }
    
    /**
     * @brief Return the `i`th segment in the path.
     *
     * If there are no complete segments in the path, return `std::nullopt`.
     */
    std::optional<BSpline<T,N>> operator[](size_t i) const {
        if (n_segments() < 1) return std::nullopt;
        return BSpline<T,N>(
            knots[i    ],
            knots[i + 1],
            knots[i + 2],
            knots[i + 3]
        );
    }
    
};


/**
 * @brief An extendable path which passes through a sequence of knots.
 * 
 * CatromSplinePaths are defined by sequence of knots, where each Catmull-Rom
 * spline segment is defined by four consecutive knots in the sequence.
 * 
 * Catmull-Rom splines are C1 continuous. The first and last knots are not
 * part of the curve.
 */
template <typename T, index_t N>
class CatromSplinePath : public SplinePath<T,N,CatromSpline<T,N>,CatromSplinePath<T,N>> {
public:
    
    using Knot = VecType<T,N>;
    
    /// The sequence of knots which the path passes through.
    std::vector<VecType<T,N>> knots;
    
    /// Construct an empty path.
    CatromSplinePath() = default;
    /// Construct a path from a sequence of knots.
    CatromSplinePath(const std::initializer_list<VecType<T,N>> knots): 
        knots(knots) {}
    /// Construct a path having a single segment.
    CatromSplinePath(const CatromSpline<T,N>& segment):
        knots {segment.control_points(), segment.control_points() + 4}
    {}
    
    /**
     * @brief Return the number of complete segments in the path.
     *
     * A single segment is defined by four consecutive control points.
     */
    size_t n_segments() const {
        size_t n = knots.size();
        if (n < 4) return 0;
        return knots.size() - 3;
    }
    
    /**
     * @brief Return the `i`th segment in the path.
     *
     * If there are no complete segments in the path, return `std::nullopt`.
     */
    std::optional<CatromSpline<T,N>> operator[](size_t i) const {
        if (n_segments() < 1) return std::nullopt;
        return CatromSpline<T,N>(
            knots[i    ],
            knots[i + 1],
            knots[i + 2],
            knots[i + 3]
        );
    }
    
};


/**
 * @brief An extendable path defined by a sequence of knots and tangent points.
 * 
 * Bezier paths are concatenated cubic Bezier splines, where the first knot
 * of each segment is the last knot of the previous segment.
 * 
 * Bezier curves are C0 continuous, and pass through the guide knots.
 * It is allowed to have "broken tangents" at each knot, though the
 * curve can be made C1 continuous by setting a tangent velocity for each knot.
 */
template <typename T, index_t N>
class BezierPath : public SplinePath<T,N,BezierSpline<T,N>,BezierPath<T,N>> {
private:
    std::vector<VecType<T,N>> knots;
    std::vector<VecType<T,N>> tangents;
    
public:
    
    /// A pair of tangents for a knot or segment in a Bezier path.
    struct Tangent {
        /// The tangent for s in the negative direction.
        VecType<T,N> t0;
        /// The tangent for s in the positive direction.
        VecType<T,N> t1;
    };
    
    /// A knot in a Bezier path, along with its tangents.
    struct Knot {
        VecType<T,N> knot;
        Tangent      tangents;
    };
    
    /// Construct an empty path.
    BezierPath() = default;
    /// Construct a path from a sequence of knots and tangents.
    BezierPath(const BezierSpline<T,N>& segment):
        knots    {segment.p0, segment.p1},
        tangents {segment.t0, segment.t1}
    {}
    
    /// Return the number of complete segments in the path.
    size_t n_segments() const {
        size_t n = knots.size();
        if (n < 2) return 0;
        return knots.size() - 1;
    }
    
    /**
     * @brief Return the `i`th segment in the path.
     *
     * If there are no complete segments in the path, return `std::nullopt`.
     */
    std::optional<BezierSpline<T,N>> operator[](size_t i) const {
        if (n_segments() < 1) return std::nullopt;
        return {
            knots   [i        ],
            tangents[2 * i    ],
            tangents[2 * i + 1],
            knots   [i     + 1]
        };
    }
    
    /// Return a reference to the `i`th knot in the path.
    VecType<T,N>& knot(size_t i)       { return knots[i]; }
    /// Return the `i`th knot in the path.
    VecType<T,N>  knot(size_t i) const { return knots[i]; }
    
    /**
     * @brief Return the tangents for the `i`th knot in the path.
     */
    Tangent tangents_for_knot(size_t i) const {
        size_t i0 = i > 0 ? (2 * i - 1) : 0;
        size_t i1 = std::min(2 * i, tangents.size());
        return {
            tangents[i0],
            tangents[i1]
        };
    }
    
    /**
     * @brief Set the tangents on either side of the `i`th knot in the path.
     *  
     * If `i` is the first or last knot, then only the inward tangent
     * is set.
     */
    void set_tangents_for_knot(size_t i, const Tangent& t) {
        if (i > 0) {
            tangents[2 * i - 1] = t.t0;
        }
        if (i < knots.size() - 1) {
            tangents[2 * i] = t.t1;
        }
    }
    
    /**
     * @brief Set the velocity for the `i`th knot in the path.
     *
     * The segments to the left and right of the knot will be made C1 continuous.
     */
    void set_velocity_for_knot(size_t i, VecType<T,N> v) {
        if (i > 0) {
            tangents[2 * i - 1] = knots[i] - v;
        }
        if (i < knots.size() - 1) {
            tangents[2 * i] = knots[i] + v;
        }
    }
    
    /**
     * @brief Return the tangents for the `i`th segment in the path.
     */
    Tangent tangents_for_segment(size_t i) const {
        return {
            tangents[2 * i],
            tangents[2 * i + 1]
        };
    }
    
    /**
     * @brief Set the tangents for the `i`th segment in the path.
     */
    void set_tangents_for_segment(size_t i, const Tangent& t) {
        tangents[2 * i    ] = t.t0;
        tangents[2 * i + 1] = t.t1;
    }
    
    /**
     * @brief Add a knot to the path, with tangents for the new segment.
     */
    void add_knot(const Knot& k) {
        knots.push_back(k.knot);
        tangents.push_back(k.t.t0);
        tangents.push_back(k.t.t1);
    }
    
    /**
     * @brief Add a knot to the path, with a symmetrical tangent for the previous knot.
     */
    void add_knot(VecType<T,N> knot, VecType<T,N> tangent) {
        size_t n = knots.size();
        VecType<T,N> t0 {(T)0};
        if (n > 0) {
            // add a symmetrical tangent
            size_t j = n - 1;
            VecType<T,N> k0 = knots[j]; 
            t0 = tangents[2 * j + 1] - k0;
        }
        knots.push_back(knot);
        tangents.push_back(t0);
        tangents.push_back(tangent);
    }
    
    /**
     * @brief Remove the `i`th knot from the path.
     *
     * If `i` is not in the range of existing knots, return `std::nullopt`.
     */
    std::optional<Knot> remove_knot(size_t i) {
        if (i >= knots.size()) return std::nullopt;
        VecType<T,N> knot = knots[i];
        Tangent t = tangents_for_knot(i);
        knots.erase(knots.begin() + i);
        tangents.erase(tangents.begin() + 2 * i, tangents.begin() + 2 * i + 2);
        return Knot {knot, t};
    }
    
    /**
     * @brief Remove all knots and their tangents from the path.
     */
    void clear() {
        knots.clear();
        tangents.clear();
    }
    
};


/**
 * @brief An extendable path defined by a sequence of knots and tangent velocities.
 *
 * Hermite paths are concatenated cubic Hermite splines, where the first knot
 * of each segment is the last knot of the previous segment. Each knot has an
 * associated tangent velocity.
 *
 * Hermite curves are C1 continuous, and pass through their knots.
 */
template <typename T, index_t N>
class HermitePath : public SplinePath<T,N,HermiteSpline<T,N>,HermitePath<T,N>> {
public:
    
    /// A knot in a Hermite path, along with its tangent velocity.
    struct Knot {
        /// The position of the knot.
        VecType<T,N> knot;
        /// The velocity of curve at the knot.
        VecType<T,N> velocity;
    };
    
    /// The sequence of knots which the path passes through.
    std::vector<Knot> knots;
    
    /// Construct an empty path.
    HermitePath() = default;
    /// Construct a path from a sequence of knots.
    HermitePath(const std::initializer_list<Knot> knots): 
        knots(knots) {}
    /// Construct a path having a single segment.
    HermitePath(const HermiteSpline<T,N>& segment):
        knots {segment.control_points(), segment.control_points() + 4}
    {}
    
    /**
     * @brief Return the number of complete segments in the path.
     *
     * Each segment is defined by four consecutive control points.
     */
    size_t n_segments() const {
        size_t n = knots.size();
        if (n < 2) return 0;
        return knots.size() - 1;
    }
    
    /**
     * @brief Return the `i`th segment in the path.
     *
     * If there are no complete segments in the path, return `std::nullopt`.
     */
    std::optional<HermiteSpline<T,N>> operator[](size_t i) const {
        if (knots.size() < 2) return std::nullopt;
        const Knot& k0 = knots[i    ];
        const Knot& k1 = knots[i + 1];
        return CatromSpline<T,N>(
            k0.knot,
            k0.velocity,
            k1.knot,
            k1.velocity
        );
    }
    
};

/// @} // addtogroup spline

} // namespace geom
