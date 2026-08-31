#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <geomc/linalg/Ray.h>
#include <geomc/linalg/Vec.h>
#include <geomc/shape/Shape.h>

// todo: remove redundant fn names.
// todo: sdf normal function
// todo: furthest point function

// Coordinates come in two flavors: positions (`point_t`, i.e. `lo`, `hi`,
// corners, containment) and differences between positions (`point_diff_t`:
// extents, displacements, distances). The latter may be negative even when the
// coordinate type cannot be, so they use `diff_t` — the signed counterpart of
// `T`. For signed and floating point `T`, `diff_t == T` and there is no
// difference; the split only matters for unsigned `T`.
//
// note: because `sdf()` and `normal()` return signed quantities, an unsigned
// Rect does not model the ProjectionObject / SdfObject concepts (whose results
// are expected to convert back to `point_t` / `elem_t`). This is deliberate:
// those functions are meaningless in unsigned space. A future cleanup may push
// this affine point/vector distinction through the rest of the library.

// all boundaries are inclusive. this is inconsistent with half-open interval convention,
// but:
// - union of rect with a point should thereafter contain the point
// - with a half-open interval, we'd need to "increment" the upper bound
//   to contain p
// - which is slow and ugly
// - the upper and lower bounds now mean different things, and
//   that means you need to increment and decrement
//   Sometimes And Not Other Times.
// - compare to: "both endpoints are symmetrical and the same logic applies everywhere"
// - the biggest downside: ownership of points on the faces of abutting Rects
//   is ambiguous (for Rects with real-valued coordinates).

namespace geom {

/****************************
 * Rect Class               *
 ****************************/

/** @ingroup shape
 *  @brief An N-dimensional axis-aligned interval.
 *
 * This class works naturally with both 1D and multi-dimensional ranges.
 * 1D ranges have point type `T`, while multi-dimensional ranges have point
 * type `Vec<T,N>`. A typedef of this may be accessed via:
 *
 *     Rect<T,N>::point_t
 *
 * Example:
 *
 * * `Rect<double,4>::point_t` is `Vec<double,4>`; `lo` and `hi` are `Vec4d`.
 * * `Rect<double,1>::point_t` is `double`; `lo` and `hi` are `double`.
 *
 * If any coordinate of `lo` is greater than the same coordinate of `hi`, then
 * the `Rect` is empty.
 *
 * Rect boundaries are inclusive. The points `lo` and `hi` are considered
 * to be inside the Rect. Therefore, `Rect`s that contain only a single point, edge,
 * or face are considered non-empty.
 *
 * Note that this convention differs slightly from conventional "interval" logic
 * on integers, wherein the upper boundary is excluded. This convention was chosen
 * so that adding any point to a Rect as a vertex ensures that the Rect thereafter
 * includes that point. To iterate over the range within an integer Rect, consider
 * using a GridIterator object, which abstracts away the boundary logic for you.
 */
template <typename T, index_t N>
class Rect: public Dimensional<T,N> {
protected:
    typedef PointType<T,N> ptype;

public:
    using typename Dimensional<T,N>::point_t;

    /**
     * @brief The signed coordinate type used for quantities produced by
     * differencing two positions (extents, displacements, and distances).
     *
     * Equal to `T` for signed and floating point coordinate types. For unsigned
     * `T`, this is the signed integer of the same width, so that a difference —
     * which may be negative even when a position cannot be — is represented
     * correctly instead of wrapping around.
     */
    using diff_t = geom::signed_type_t<T>;
    /// The point type for a difference between two positions; a `Vec<diff_t,N>`
    /// if `N > 1`, otherwise a `diff_t`.
    using point_diff_t = VecType<diff_t,N>;

protected:
    // point type helper for difference-space quantities
    typedef PointType<diff_t,N> dtype;

    // Promote a position into the signed "difference" space, so that
    // subtractions do not wrap for unsigned coordinate types. This is a no-op
    // (identity) whenever `diff_t == T`.
    static point_diff_t _to_diff(point_t p) {
        return (point_diff_t) p;
    }

    // Demote a difference-space quantity back to a position. Callers are
    // responsible for ensuring the value is representable (non-negative for
    // unsigned `T`); out-of-range results wrap, as with any narrowing.
    static point_t _from_diff(point_diff_t p) {
        return (point_t) p;
    }

public:

    /// Lower extremes
    point_t lo;
    /// Upper extremes
    point_t hi;

    static const point_t endpoint_measure;

public:


    /**
     * @brief Construct an empty interval.
     *
     * Sets lower and upper bounds to the maximum and minimum values of `T`, respectively.
     *
     * A union between this Rect and any finite shape is an identity operation.
     */
    constexpr Rect():
        lo(std::numeric_limits<T>::max()),
        hi(std::numeric_limits<T>::lowest()) {}

    /**
     * @brief Construct a Rect with extremes `lo` and `hi`.
     * If for any axis `lo > hi`, the Rect is empty.
     *
     * @param lo Lower extreme
     * @param hi Upper extreme
     */
    constexpr Rect(point_t lo, point_t hi):
        lo(lo),
        hi(hi) {}

    /**
     * @brief Construct a Rect containing only the point `p`.
     */
    explicit constexpr Rect(point_t p):
        lo(p),
        hi(p) {}
    
    /**
     * @brief Construct a Rect from lower-dimensional Rects.
     */
    template <index_t... J>
    requires ((J + ...) == N)
    constexpr Rect(const Rect<T,J>&... r):
        lo(r.lo...),
        hi(r.hi...) {}


    /****************************
     * Static members           *
     ****************************/

    /// A Rect covering the range [0, 1] along all axes.
    static const Rect<T,N> unit_interval;

    /// A Rect covering the range [-1, 1] along all axes.
    static const Rect<T,N> signed_unit_interval;

    /// A Rect that contains all points.
    static const Rect<T,N> full;

    /// A Rect that contains no points.
    static const Rect<T,N> empty;

    /**
     * @brief Construct a Rect from a center point and extent.
     * @param c Center of the new Rect
     * @param dims Lengths of each axis; may be negative.
     * @return A new Rect with its center at `c`.
     */
    inline static Rect<T,N> from_center(point_t c, point_diff_t dims) {
        point_diff_t c_d = _to_diff(c);
        point_diff_t em  = _to_diff(endpoint_measure);
        return Rect<T,N>::from_corners(
            _from_diff(c_d -  dims       / 2),
            _from_diff(c_d + (dims - em) / 2)
        );
    }

    /**
     * @brief Construct a Rect from a corner point and an extent.
     * @param edge Arbitrary corner point.
     * @param dims Lengths of each axis relative to given corner. Lengths may be negative.
     * @return A new Rect with one corner at `corner`.
     */
    inline static Rect<T,N> from_edge(point_t edge, point_diff_t dims) {
        point_diff_t em = _to_diff(endpoint_measure);
        return Rect<T,N>::from_corners(
            edge,
            _from_diff(_to_diff(edge) + dims - em)
        );
    }

    /**
     * @brief Construct a Rect containing the two corners `c1` and `c2`.
     * @param c1 A corner of the Rect
     * @param c2 The corner opposite `c1`
     * @return A new Rect with corners at `c1` and `c2`.
     */
    inline static Rect<T,N> from_corners(
            point_t c1,
            point_t c2)
    {
        return Rect<T,N>(
            std::min(c1, c2),
            std::max(c1, c2)
        );
    }

    /**
     * @brief Construct a Rect at the origin having dimensions `dims`.
     *
     * @param dims Lengths of each axis.
     * @return Rect<T,N> A new Rect with its center at the origin.
     */
    inline static Rect<T,N> from_size(point_diff_t dims) {
        return Rect<T,N>::from_center(point_t {}, dims);
    }

    /**
     * @brief Construct a Rect containing all the points in the sequence
     * between `begin` and `end`.
     * @tparam PointIterator An iterator to a `point_t`.
     * @param begin Iterator to the first point.
     * @param end Iterator to the off-end point.
     * @return A new Rect containing all the points in the given sequence.
     */
    template <typename PointIterator>
    inline static Rect<T,N> from_point_sequence(PointIterator begin, PointIterator end) {
        Rect<T,N> r;
        for (PointIterator i = begin; i != end; ++i) {
            r |= *i;
        }
        return r;
    }

    /**
     * @brief Test whether a point is in the N-dimensional range `[lo, hi]`.
     *
     * If `lo` > `hi` along any axis, then the range is empty and the
     * function returns `false`.
     *
     * @param lo Lower extreme
     * @param hi Upper extreme
     * @param pt Test point
     * @return `true` if `pt` is inside the exremes
     */
    static bool contains(VecType<T,N> lo, VecType<T,N> hi, VecType<T,N> pt) {
        for (index_t axis = 0; axis < N; axis++) {
            T v = coord(pt, axis);
            if (v < coord(lo, axis) or v > coord(hi, axis)) {
                return false;
            }
        }
        return true;
    }
    
    static constexpr bool admits_cusps() { return false; }

    /*****************************
     * Operators                 *
     *****************************/

    /**
     * @brief Dimensionally slice a Rect.
     *
     * Return a Rect of dimension `M` by taking the range for each axis given by
     * the arugments. For example, the result of `r(0,2)` will be a 2D Rect
     * with the x and z ranges of `this` Rect.
     */
    template <typename... I>
    requires (... and std::convertible_to<I, size_t>)
#if __cplusplus >= 202302L
    [[deprecated("use operator[] instead")]]
#endif
    Rect<T,sizeof...(I)> operator()(I... indices) const {
        constexpr index_t M = sizeof...(I);
        Rect<T,M> out {
            point_t{ coord(lo, indices)...},
            point_t{ coord(hi, indices)...}
        };
        return out;
    }

#if __cplusplus >= 202302L
    /**
     * @brief Dimensionally slice a Rect.
     *
     * Return a Rect of dimension `M` by taking the range for each axis given by
     * the arugments. For example, the result of `r[0,2]` will be a 2D Rect
     * with the x and z ranges of `this` Rect.
     */
    template <typename... I>
    requires (... and std::convertible_to<I, size_t>)
    Rect<T,sizeof...(I)> operator[](I... indices) const {
        constexpr index_t M = sizeof...(I);
        Rect<T,M> out {
            point_t{ coord(lo, indices)...},
            point_t{ coord(hi, indices)...}
        };
        return out;
    }
#endif

    /**
     * @brief Interval union.
     *
     * @return A Rect fully containing `this` and box `b`.
     */
    Rect<T,N> operator|(const Rect<T,N> &b) const {
        return Rect<T,N>(std::min(lo, b.lo),
                         std::max(hi, b.hi));
    }

    /**
     * @brief Interval union.
     *
     * Extend `this` to fully contain `b`.
     * @return A reference to `this`, for convenience.
     */
    Rect<T,N>& operator|=(const Rect<T,N> &b) {
        // box union
        hi = std::max(b.hi, hi);
        lo = std::min(b.lo, lo);
        return *this;
    }

    /**
     * @brief Point union.
     *
     * @return A Rect fully containing `this` and the point `p`.
     */
    Rect<T,N> operator|(const point_t& p) {
        // point union
        return Rect<T,N>(
            std::min(lo, p),
            std::max(hi, p)
        );
    }

    /**
     * @brief Point union.
     *
     * Extend `this` to fully contain `p`.
     * @return a reference to `this`, for convenience.
     */
    Rect<T,N>& operator|=(const point_t& p) {
        // point union
        lo = std::min(lo, p);
        hi = std::max(hi, p);
        return *this;
    }

    /**
     * @brief Interval intersection.
     *
     * Compute a Rect representing the area overlapped by both `this` and `b`.
     * This is the unsigned intersection, which means that an intersection with an
     * empty interval yields another empty interval. If empty intervals should
     * instead subtract from nonempty intervals, use signed_intersection().
     *
     * The unsigned intersection is the most common form of interval intersection,
     * and is the more efficient of the two.
     */
    Rect<T,N> operator&(const Rect<T,N> &b) const {
        return Rect<T,N>(std::max(lo, b.lo),
                         std::min(hi, b.hi));
    }

    /**
     * @brief Interval intersection.
     *
     * Confine this Rect to the area overlapping `b`. This is the unsigned
     * intersection, which means that an intersection with an empty interval yields
     * another empty interval. If empty intervals should instead subtract from
     * nonempty intervals, use signed_intersection().
     *
     * The unsigned intersection is the most common form of interval intersection,
     * and is the more efficient of the two.
     * @return A reference to `this` for convenience.
     */
    Rect<T,N>& operator&=(const Rect<T,N> &b) {
        lo = std::max(lo, b.lo);
        hi = std::min(hi, b.hi);
        return *this;
    }

    /**
     * @brief Translation.
     *
     * @param dx Amount by which to translate this region.
     *   Add `dx` to the coordinates of all the bounds.
     * @return A translated Rect.
     */
    Rect<T,N> operator+(point_diff_t dx) const {
        return Rect<T,N>(
            _from_diff(_to_diff(lo) + dx),
            _from_diff(_to_diff(hi) + dx)
        );
    }

    /**
     * @brief Translation.
     *
     * @param dx Amount by which to translate this region.
     *   Subtract `dx` from the coordinates of all the bounds.
     * @return A translated Rect.
     */
    Rect<T,N> operator-(point_diff_t dx) const {
        return Rect<T,N>(
            _from_diff(_to_diff(lo) - dx),
            _from_diff(_to_diff(hi) - dx)
        );
    }


    /**
     * @brief Translation.
     *
     * @param dx Amount by which to translate this region.
     *   Add `dx` to the coordinates of all the bounds.
     * @return  A reference to `this`, for convenience.
     */
    Rect<T,N>& operator+=(point_diff_t dx) {
        lo = _from_diff(_to_diff(lo) + dx);
        hi = _from_diff(_to_diff(hi) + dx);
        return *this;
    }


    /**
     * @brief Translation.
     *
     * @param dx Amount by which to translate this region.
     *   Subtract `dx` from the coordinates of all the bounds.
     * @return  A reference to `this`, for convenience.
     */
    Rect<T,N>& operator-=(point_diff_t dx) {
        lo = _from_diff(_to_diff(lo) - dx);
        hi = _from_diff(_to_diff(hi) - dx);
        return *this;
    }


    /**
     * @brief Interval exclusion.
     *
     * For any dimension where `this` fully contains `other`, the result is the
     * empty interval. Any dimension of `other` that is empty has no effect.
     *
     * @param other Range to be excluded from the range of this `Rect`.
     * @return A new `Rect` whose range does not overlap `other`.
     */
    Rect<T,N> operator-(const Rect<T,N>& other) {
        // todo: instead of returning empty() when an interval is split in two,
        //   choose the larger
        // todo: make an exclude_signed() and handle empty-interval case,
        //   and make it behave sensibly? e.g.:
        //   •   a  - (-b) =   a + b
        //   • (-a) -   b  = -(a + b)
        //   • (-a) - (-b) = -(a - b)
        Rect<T,N> out;
        for (index_t i = 0; i < N; ++i) {
            T b_lo = coord(other.lo, i);
            T a_lo = coord(lo, i);
            T a_hi = coord(hi, i);
            T b_hi = coord(other.hi, i);
            if (b_lo > b_hi) continue;
            bool pick_lo = b_lo < a_lo;
            bool pick_hi = b_hi > a_hi;
            coord(lo, i) = pick_lo ? a_hi : std::min(a_hi, b_lo);
            coord(hi, i) = pick_hi ? a_lo : std::max(b_hi, a_lo);
        }
        return out;
    }


    /**
     * @brief Interval exclusion.
     *
     * For any dimension where `this` fully contains `other`, the result is the
     * empty interval. Any dimension of `other` that is empty has no effect.
     *
     * @param other Range to remove from this `Rect`.
     * @return A reference to `this`, for convenience.
     */
    Rect<T,N>& operator-=(const Rect<T,N>& other) {
        *this = (*this) - other;
        return *this;
    }


    /**
     * @brief Equality test.
     *
     * @return `true` if and only if all the corresponding extremes of `b` are the same.
     */
    bool operator==(const Rect<T,N>& b) const {
        return (hi == b.hi) and (lo == b.lo);
    }

    /**
     * @brief Inequality test.
     *
     * @return `true` if and only if any extreme of `b` is different from
     * the corresponding extreme in `this`.
     */
    bool operator!=(const Rect<T,N>& b) const {
        return (hi != b.hi) or (lo != b.lo);
    }

    /**
     * @brief Scale transformation.
     *
     * @param a Scale factor.
     * @return A new Rect, scaled about the origin by factor `a`.
     */
    Rect<T,N> operator*(point_t a) const {
        // a sign change can swap lo and hi for any coordinate.
        // we have to do max and min again after the multiply:
        point_t lo_a = lo * a;
        point_t hi_a = hi * a;
        return Rect<T,N>(
            std::min(lo_a, hi_a),
            std::max(lo_a, hi_a));
    }

    /**
     * @brief Scale transformation.
     *
     * Scale this Rect about the origin by factor `a`.
     *
     * @param a Scale factor
     * @return A reference to `this`, for convenience.
     */
    Rect<T,N>& operator*=(point_t a) {
        point_t lo_a = lo * a;
        point_t hi_a = hi * a;
        lo = std::min(lo_a, hi_a);
        hi = std::max(lo_a, hi_a);
        return *this;
    }

    /**
     * @brief Scale transformation.
     *
     * @param a Scale factor.
     * @return A new Rect, scaled about the origin by multiple `1 / a`.
     */
    Rect<T,N> operator/(point_t a) const {
        auto lo_a = lo / a;
        auto hi_a = hi / a;
        return {
            std::min(lo_a, hi_a),
            std::max(lo_a, hi_a)
        };
    }

    /**
     * @brief Scale transformation.
     *
     * Scale this Rect about the origin by factor `1 / a`.
     * @param a Scale factor
     * @return A reference to `this`, for convenience.
     */
    Rect<T,N>& operator/=(point_t a) {
        *this = (*this) / a;
        return *this;
    }
    
    /**
     * @brief Negate all coordinates.
     *
     * The volume of the box keeps its sign; i.e., if it is nonempty,
     * its negation is nonempty.
     */
    Rect<T,N> operator-() const {
        return {-hi, -lo};
    }

    /**
     * @brief Interval cartesian product.
     *
     * Extrude this Rect into a higher dimension, by treating the extents of
     * `r` as the extents along the new dimensions. In other words, concatenate
     * the coordinates of `this` and `r` into a new Rect.
     *
     * For example, the product of a 2D rectangle with a 1D range is a 3D box.
     * If multiplied in that order, the product will have the extents of the rectangle
     * along its first two axes and the extents of the 1D range along its third axis.
     *
     * @tparam M Dimensionality of `r`.
     */
    template <index_t M>
    Rect<T,M+N> operator*(const Rect<T,M>& r) const {
        typedef PointType<T,M> r_ptype;
        Rect<T,M+N> o;

        std::copy(  ptype::iterator(lo),     ptype::iterator(lo)   + N, o.lo.begin());
        std::copy(  ptype::iterator(hi),     ptype::iterator(hi)   + N, o.hi.begin());
        std::copy(r_ptype::iterator(r.lo), r_ptype::iterator(r.lo) + M, o.lo.begin() + N);
        std::copy(r_ptype::iterator(r.hi), r_ptype::iterator(r.hi) + M, o.hi.begin() + N);

        return o;
    }

    /**
     * @brief Element-wise typecast.
     *
     * @return A new Rect, with elements all of type `U`.
     */
    template <typename U, index_t M>
    operator Rect<U,M>() const {
        if constexpr (M == N) {
            return Rect<U,M>(
                VecType<U,N>(lo),
                VecType<U,N>(hi)
            );
        } else {
            return Rect<U,M>(
                VecType<U,N>(lo).template resized<M>(),
                VecType<U,N>(hi).template resized<M>()
            );
        }
    }


    /*****************************
     * Public Methods            *
     *****************************/

    /**
     * @brief Point containment test.
     * @return `true` if and only if `pt` is inside this rectangle.
     * Points on the surface of the Rect are considered to be contained by it.
     */
    bool contains(point_t pt) const {
        for (index_t axis = 0; axis < N; axis++) {
            T v = coord(pt, axis);
            if (v < coord(lo, axis) or v > coord(hi, axis)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Box containment.
     *
     * @return `true` if and only if this `Rect` *fully contains* `box`, in other
     * words that there is no point contained by `box` which is not contained by
     * `this`.
     */
    bool contains(const Rect<T,N>& box) const {
        for (index_t axis = 0; axis < N; axis++) {
            // inner box must not protrude along any axis
            if (coord(box.lo, axis) < coord(lo, axis) or
                coord(box.hi, axis) > coord(hi, axis))
            {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Return the one dimensional range spanned along axis `k`.
     */
    Rect<T,1> axis(size_t k) const {
        return {
            coord(lo, k),
            coord(hi, k)
        };
    }

    /**
     * @brief Range intersection test.
     * @return `true` if and only if there is a point overlapped by both `Rect`s.
     */
    bool intersects(const Rect<T,N>& box) const {
        for (index_t axis = 0; axis < N; axis++) {
            // disjoint on this axis?
            if (coord(hi, axis)     <= coord(box.lo, axis) or
                coord(box.hi, axis) <= coord(lo, axis)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Transformed range intersection test.
     *
     * Alias for `box.intersects(*this)`.
     */
    bool intersects(const Transformed<Rect<T,N>>& box) const {
        // delegate to OrientedRect
        return box.intersects(*this);
    }

    bool intersects(const Sphere<T,N>& sph) const {
        return sph.intersects(*this);
    }
    
    template <ConvexObject Shape>
    requires (Shape::N == N) and std::same_as<T, typename Shape::elem_t>
    bool intersects(const Shape& other) const {
        return geom::intersects(
            as_any_convex(*this),
            as_any_convex(other)
        );
    }

    /**
     * @brief Center point.
     * @return The center point of this region.
     */
    point_t center() const {
        // computed as lo + (hi - lo) / 2 rather than (lo + hi) / 2: the sum
        // (lo + hi) overflows near the top of the type's range and is not
        // representable at all for unsigned coordinates whose sum exceeds the
        // maximum. The lo + (hi - lo)/2 form is overflow-safe for a non-empty
        // Rect (hi >= lo, so hi - lo is representable) in the coordinate type
        // itself, without needing a wider or signed type.
        return lo + (hi - lo + endpoint_measure) / 2;
    }

    /**
     * @brief Axial size.
     * @return The size of this region along each axis.
     *
     * Note that for integer type Rects, since both the high and low boundaries
     * are included, the length along each axis is `hi - lo + 1`.
     *
     * The result is a difference of positions, and so has signed type
     * `point_diff_t`. An empty Rect (with `lo > hi` on some axis) reports a
     * negative extent along that axis.
     */
    point_diff_t dimensions() const {
        return _to_diff(hi) - _to_diff(lo) + _to_diff(endpoint_measure);
    }

    /**
     * @brief Change the size of the region, adjusting about its center.
     *
     * A negative length is treated as its absolute value.
     *
     * @param dim New lengths along each axis.
     */
    void set_dimensions(point_diff_t dim) {
        // rule: preserve the center point.
        //       for odd-length boxes, this is the center element.
        //       for even-length boxes, it is the middle zero-indexed "fencepost".
        // [0|1|2|3]
        //   [1|2|3]

        // [0|1|2|3]
        //   [1|2]

        // [0|1|2|3|4]
        //   [1|2|3]

        // [0|1|2|3|4]
        // [0|1|2|3]

        // normalize the sign per-axis (Vec::abs() is unavailable for the N == 1
        // scalar case, so clamp element-wise):
        for (index_t i = 0; i < N; ++i) {
            coord(dim, i) = std::abs(coord(dim, i));
        }
        // total change to distribute across the two sides. split so that the
        // sum applied to lo and hi is exactly `delta` even for odd deltas
        // (integer division would otherwise lose the odd unit); the extra unit
        // goes to the hi side, consistent with the "hi fencepost" center rule.
        point_diff_t delta   = dim - dimensions();
        point_diff_t lo_grow = delta / 2;
        point_diff_t hi_grow = delta - lo_grow;
        lo = _from_diff(_to_diff(lo) - lo_grow);
        hi = _from_diff(_to_diff(hi) + hi_grow);
    }
    
    /**
     * @brief Return the `i`th corner point of the Rect.
     */
    point_t corner(index_t i) const {
        point_t c;
        for (index_t j = 0; j < N; ++j) {
            coord(c, j) = (i & (1 << j)) ? coord(hi, j) : coord(lo, j);
        }
        return c;
    }

    /**
     * @brief Reconstruct from corner points.
     *
     * Re-configure this region to exactly contain the two given points.
     */
    void set_corners(point_t corner1, point_t corner2) {
        hi = std::max(corner1, corner2);
        lo = std::min(corner1, corner2);
    }

    /**
     * @brief Place a copy at a new location.
     *
     * Make a new Rect having the same size as this one, but with its center
     * at `c`.
     * @param c New center point.
     */
    Rect<T,N> centered_on(point_t c) const {
        return Rect<T,N>::from_center(c, dimensions());
    }

    /**
     * @brief Normalize the sign of the region.
     *
     * If any extent has negative sign (lo > hi), the `lo` and `hi` coordinates
     * for that axis are swapped. The result is a non-empty Rect.
     */
    Rect<T,N> abs() const {
        Rect<T,N> out;
        for (index_t i = 0; i < N; ++i) {
            coord(out.lo, i) = std::min(coord(lo, i), coord(hi, i));
            coord(out.hi, i) = std::max(coord(lo, i), coord(hi, i));
        }
        return out;
    }

    /**
     * @brief Return a rect having the same aspect (ratios of axial lengths) as this one,
     * sized and positioned to minimally contain `other`.
     *
     * The fitted box will be a rigid rescaling and translation of this box. No reflections
     * or rotations will be performed. The signs of the dimensions will be the same as this
     * Rect.
     *
     * Not available for integral ranges, since in general the aspect can't be preserved.
     *
     * Zero dimensions will be ignored and not adjusted. To expand the zero dimensions, take
     * the union of `other` with the result.
     *
     * @param other Rect to contain.
     */
    Rect<T,N> fitted_to(Rect<T,N> other) const requires (not std::integral<T>) {
        T max_ratio = 0;
        point_t d0 = dimensions();
        point_t d1 = other.dimensions();
        // find the dimension which constrains the final size:
        for (index_t i = 0; i < N; ++i) {
            // do not try to "expand" any dimensions which are empty:
            auto ratio = (coord(d0, i) == 0) ? 0 : (coord(d1, i) / coord(d0, i));
            max_ratio = std::max(max_ratio, std::abs(ratio));
        }
        return Rect<T,N>::from_center(other.center(), d0 * max_ratio);
    }

    /**
     * @brief Morphological dilation.
     *
     * Return a Rect with the boundary extended coordinate-wise by the amount `c`
     * in all directions. A negative `c` erodes (shrinks) the Rect.
     */
    Rect<T,N> dilated(point_diff_t c) const {
        return Rect<T,N>(
            _from_diff(_to_diff(lo) - c),
            _from_diff(_to_diff(hi) + c)
        );
    }

    /**
     * @brief Minkowski sum.
     *
     * Combine this Rect with another using the Minkowski sum.
     */
    Rect<T,N> minkowski_sum(const Rect<T,N>& other) const {
        return Rect<T,N>(lo + other.lo, hi + other.hi);
    }

    /**
     * @brief N-dimensional volume.
     *
     * Compute the measure of the region contained by this Rect as the product
     * of its dimensions. Result measures a length if `N` is 1; an area if `N` is 2; a
     * volume if 3; etc.
     *
     * The volume of an empty or negative Rect is 0.
     *
     * @return The volume of this region.
     */
    T measure_interior() const {
        // accumulate in signed space so the negative-extent clamp works for
        // unsigned `T`; the clamped product is non-negative and converts back
        // to `T` losslessly.
        diff_t vol = 1;
        point_diff_t dim = dimensions();
        for (index_t axis = 0; axis < N; axis++) {
            vol *= std::max(coord(dim, axis), (diff_t) 0);
        }
        return (T) vol;
    }
    
    /**
     * @brief Measure of the boundary of this Rect.
     * 
     * Measures the length/perimter/surface area of this Rect.
     */
    T measure_boundary() const {
        return 2 * face_areas().sum();
    }
    
    /**
     * @brief Measure of the axial faces of this Rect.
     *
     * The `i`th component of the result is the area of the face perpendicular to the `i`th axis.
     */
    point_diff_t face_areas() const {
        point_diff_t axial_areas;
        diff_t prod = 1;
        point_diff_t dims = dimensions();
        // include all the extents preceding this axis (but not the axis itself)
        for (index_t i = 0; i < N - 1; ++i) {
            coord(axial_areas, i) = prod;
            prod *= coord(dims, i);
        }
        // include all the extents following this axis (but not the axis itself)
        prod = 1;
        for (index_t i = N - 1; i >= 0; --i) {
            coord(axial_areas, i) *= prod;
            prod *= coord(dims, i);
        }
        return axial_areas;
    }
    
    /**
     * @brief Return a "normalized" length measure encoding the size of this box.
     *
     * Computed as the `N`th root of the volume.
     */
    auto length() const {
        return std::pow(measure_interior(), 1/(double)N);
    }
    
    /**
     * @brief Empty region test.
     * @return `true` if and only if this region contains no points; which
     * will be the case if `lo[i] > hi[i]` for any axis `i`.
     */
    bool is_empty() const {
        for (index_t axis = 0; axis < N; axis++) {
            T hi_i = coord(hi, axis);
            T lo_i = coord(lo, axis);
            if (hi_i < lo_i) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Clamp the coordinates of `p` to lie within this Rect.
     *
     * Result can be considered the point nearest to `p` contained in this `Rect`.
     */
    point_t clip(point_t p) const {
        return std::min(hi, std::max(lo, p));
    }

    /**
     * @brief Return the point on the surface of the shape which is nearest to `p`.
     *
     * Distinct from clip() in that `p` is projected to the boundary of
     * the Rect, regardless of whether it's inside or outside. By contrast,
     * `clip()` leaves `p` unchanged if it lies in the Rect's interior.
     */
    point_t project(point_t p) const {
        if (contains(p)) {
            // project to the closest face;
            // i.e. find the axis along which p is closest to a boundary and 
            // snap that coordinate to the boundary.
            index_t best_axis    = 0;
            T      winning_coord = 0;
            // distances are differences of positions; compute them signed so
            // they do not wrap for unsigned `T`.
            diff_t best_distance = std::numeric_limits<diff_t>::max();
            for (index_t i = 0; i < N; ++i) {
                T   p_i = coord(p,  i);
                T  lo_i = coord(lo, i);
                T  hi_i = coord(hi, i);
                diff_t to_lo = std::abs((diff_t) p_i - (diff_t) lo_i);
                diff_t to_hi = std::abs((diff_t) p_i - (diff_t) hi_i);
                diff_t  dist = std::min(to_lo, to_hi);
                if (dist < best_distance) {
                    winning_coord = to_lo < to_hi ? lo_i : hi_i;
                    best_axis     = i;
                    best_distance = dist;
                }
            }
            coord(p, best_axis) = winning_coord;
            return p;
        } else {
            return clip(p);
        }
    }
    
    /**
     * @brief Outward-facing direction.
     *
     * Returns a signed direction vector (`point_diff_t`), which for unsigned `T`
     * is not representable as a `point_t`.
     */
    point_diff_t normal(point_t p) const {
        // compute the outward-pointing normal from `p`.
        // we don't just do `p - project(p)`, because that behaves degenerately near
        // the surface of the range (where we are most likely to care about the normal).

        // nearest point in the range to `p`
        point_t clipped;
        // how many axes are outside the bound?
        index_t outside = 0;
        // which axis' plane is closest to p?
        index_t nearest_axis;
        // what is the distance to the nearest bounding plane? (a difference of
        // positions, so signed to avoid wrapping for unsigned `T`)
        diff_t nearest_plane = std::numeric_limits<diff_t>::max();
        // is the nearest plane a lower (-1) or upper (+1) bound?
        index_t nearest_sign;
        // the last axis we encountered that we were not contained by
        index_t exterior_axis;
        // is the nearest exterior face a lower bound (-1) or upper bound (+1)?
        index_t exterior_sign;

        for (index_t i = 0; i < N; ++i) {
            T x    = coord(p,  i);
            T lo_i = coord(lo, i);
            T hi_i = coord(hi, i);
            coord(clipped, i) = std::min(hi_i, std::max(lo_i, x));
            bool is_out = x < lo_i or x > hi_i;
            outside += is_out;
            // compute distances to the two bounding planes
            diff_t d0 = std::abs((diff_t) x - (diff_t) lo_i);
            diff_t d1 = std::abs((diff_t) x - (diff_t) hi_i);
            diff_t d_min = std::min(d0, d1);
            if (d_min < nearest_plane) {
                nearest_axis  = i;
                nearest_plane = d_min;
                nearest_sign  = d0 < d1 ? -1 : 1;
            }
            if (is_out) {
                exterior_axis = i;
                exterior_sign = d0 < d1 ? -1 : 1;
            }
        }
        if (outside == 0) {
            // we're contained by the range. return a normal pointing toward the nearest face.
            point_diff_t n;
            coord(n, nearest_axis) = nearest_sign;
            return n;
        } else if (outside == 1) {
            // we're outside the range, but project to an axial face. return its normal
            point_diff_t n;
            coord(n, exterior_axis) = exterior_sign;
            return n;
        } else {
            // the nearest point is on a corner / edge.
            // draw a vector from that point to `p` and unitize.
            // this still suffers numerical instability near the surface, but it's
            // near a degenerate feature, so there's not really anything we can do.
            return (_to_diff(p) - _to_diff(clipped)).unit();
        }
    }

    /**
     * @brief Return the signed distance to the surface of the shape.
     *
     * Negative inside the Rect; positive outside. The result is signed
     * (`diff_t`) and, for unsigned `T`, is not representable as an `elem_t`.
     */
    diff_t sdf(point_t p) const {
        diff_t sign = contains(p) ? -1 : 1;
        // compute the distance in signed space so it does not wrap for unsigned `T`.
        point_diff_t d = _to_diff(project(p)) - _to_diff(p);
        return sign * dtype::mag(d);
    }

    /**
     * @brief Squared distance to the interior of the Rect.
     * @return The square of the distance to the nearest point contained by this `Rect`;
     * zero if `p` is inside.
     */
    diff_t dist2(point_t p) const {
        // the displacement to the nearest contained point is a difference of
        // positions; compute it signed so it does not wrap for unsigned `T`.
        point_diff_t d = _to_diff(p) - _to_diff(clip(p));
        return dtype::mag2(d);
    }

    /**
     * @brief Map the unit interval to the region.
     *
     * Remap `s` on the `[0,1]` interval to the extents of this Rect.
     * In other words, interpolate the corners of this Rect using s as an interpolation
     * parameter.
     *
     * Inverse operation of `unmap()`.
     *
     * Values of `s` between 0 and 1 correspond to points inside this Rect.
     *
     * @see remap_transform()
     */
    point_t remap(point_t s) const {
        // todo: template for integer type? beware endpoint measure.
       return (point_t(1) - s) * lo + s * hi;
    }

    /**
     * @brief Find `p`'s fractional position within this Rect.
     *
     * Inverse operation of `remap()`.
     *
     * The numerator and denominator are differences of positions, so the result
     * is computed in signed space and returned as `point_diff_t`.
     *
     * @see unmap_transform()
     */
    point_diff_t unmap(point_t p) const {
        return (_to_diff(p) - _to_diff(lo)) / (_to_diff(hi) - _to_diff(lo));
    }
    
    /**
     * @brief Map the unit interval to the region, broadcasting to a different dimension.
     *
     * Requires that `N == 1`.
     *
     * Each coordinate of `s` is treated as a fraction within this 1-D Rect.
     */
    template <index_t M>
    requires (N == 1 and M > 1)
    Vec<T,M> remap(Vec<T,M> s) const {
        return (point_t(1) - s) * lo + s * hi;
    }
    
    /**
     * @brief Find `p`'s fractional position within this Rect, broadcasting to a different dimension.
     *
     * Requires that `N == 1`.
     *
     * Each coordinate of the result is the fraction of the corresponding coordinate of `p`
     * within this 1-D Rect.
     */
    template <index_t M>
    requires (N == 1 and M > 1)
    Vec<diff_t,M> unmap(Vec<T,M> p) const {
        diff_t lo_d = (diff_t) lo;
        return ((Vec<diff_t,M>) p - lo_d) / ((diff_t) hi - lo_d);
    }

    /**
     * @brief Return a vector which moves `this` into precise disjoint contact
     * with `other`, and a boolean indicating whether the two shapes currently
     * overlap.
     */
    std::pair<point_diff_t,bool> contact_vector(const Rect<T,N>& other) const {
        // the contact vector is a displacement, so its components are signed
        // (`diff_t`); they may be negative even for unsigned `T`.
        point_diff_t v;
        bool overlap[N];
        bool all_overlap = true;
        index_t shortest_axis = 0;
        diff_t  shortest_val  = std::numeric_limits<diff_t>::max();
        for (index_t i = 0; i < N; ++i) {
            Rect<T,1> r0 = this->axis(i);
            Rect<T,1> r1 = other.axis(i);
            diff_t d0 = (diff_t) r1.lo - (diff_t) r0.hi;
            diff_t d1 = (diff_t) r1.hi - (diff_t) r0.lo;
            overlap[i] = r0.intersects(r1);
            all_overlap = all_overlap and overlap[i];
            diff_t s = std::abs(d0) < std::abs(d1) ? d0 : d1;
            coord(v, i) = s;
            if (std::abs(s) < std::abs(shortest_val)) {
                shortest_val  = s;
                shortest_axis = i;
            }
        }
        if (all_overlap) {
            // if the two rects overlap, move along the axis that is closest to the
            // surface. once those ranges are disjoint, it doesn't matter that the
            // other axes overlap. (and separating other axes will just increase the
            // displacement distance)
            v = {};
            coord(v, shortest_axis) = shortest_val;
        } else {
            // if the two rects don't overlap, then close the gap on all
            // non-overlapping axes
            for (index_t i = 0; i < N; ++i) {
                if (overlap[i]) {
                    // overlap on this axis.
                    coord(v, i) = 0;
                }
            }
        }
        return {v, all_overlap};
    }

    point_t convex_support(point_t d) const {
        point_t o;
        for (index_t i = 0; i < N; i++) {
            T a = coord(d, i);
            coord(o, i) = coord(a < 0 ? lo : hi, i);
        }
        return o;
    }

    /// Ray-shape intersection test.
    Rect<T,1> intersect(const Ray<T,N>& r) const {
        Rect<T,1> interval = Rect<T,1>::full;
        for (index_t axis = 0; axis < N; ++axis) {
            T dx   = coord(r.direction, axis);
            T o_i  = coord(r.origin,    axis);
            T lo_i = coord(lo, axis);
            T hi_i = coord(hi, axis);
            if (dx == 0) {
                // ray direction is tangent to the test planes; no intersection on this axis.
                if (o_i < lo_i or o_i > hi_i) {
                    // origin outside the test slab; miss
                    return Rect<T,1>();
                }
            } else {
                // the ray parameter is a difference of positions over the
                // direction; compute it in signed space so it does not wrap for
                // unsigned `T`. (The parametric interval is nonetheless stored
                // as `T`; ray casts against unsigned Rects remain degenerate.)
                diff_t dx_d = (diff_t) dx;
                interval &= Rect<T,1>::from_corners(
                    ((diff_t) hi_i - (diff_t) o_i) / dx_d,
                    ((diff_t) lo_i - (diff_t) o_i) / dx_d);
            }
        }
        return interval;
    }

    Rect<T,N> bounds() const {
        return *this;
    }

}; // end Rect class


namespace detail {

template <typename T>
constexpr T lowest() {
    if constexpr (std::numeric_limits<T>::has_infinity) {
        return -std::numeric_limits<T>::infinity();
    } else {
        return std::numeric_limits<T>::lowest();
    }
}

template <typename T>
constexpr T max() {
    if constexpr (std::numeric_limits<T>::has_infinity) {
        return std::numeric_limits<T>::infinity();
    } else {
        return std::numeric_limits<T>::max();
    }
}

} // end namespace detail


// `T` and `N` are deduced from the Rect; the displacement's type follows from
// `T` in a non-deduced context so that it accepts the signed difference type.
template <typename T, index_t N>
inline Rect<T,N> operator+(VecType<signed_type_t<T>,N> dx, const Rect<T,N>& r) {
    return r + dx;
}

template <typename T, index_t N>
inline Rect<T,N> operator-(VecType<signed_type_t<T>,N> dx, const Rect<T,N>& r) {
    return r - dx;
}


// these must be declared outside the class body because
// the class definition must be complete in order to invoke
// the constexpr constructor:

template <typename T, index_t N>
constexpr Rect<T,N> Rect<T,N>::unit_interval = Rect<T,N>((T)0, (T)1);

template <typename T, index_t N>
constexpr Rect<T,N> Rect<T,N>::signed_unit_interval = Rect<T,N>((T)-1, (T)1);

template <typename T, index_t N>
constexpr Rect<T,N> Rect<T,N>::full = Rect<T,N>(
    detail::lowest<T>(),
    detail::max<T>()
);

template <typename T, index_t N>
constexpr Rect<T,N> Rect<T,N>::empty = Rect<T,N>(
    detail::max<T>(),
    detail::lowest<T>()
);

template <typename T, index_t N>
constexpr typename Rect<T,N>::point_t Rect<T,N>::endpoint_measure =
    std::is_integral<T>::value ? 1 : 0;


// hashing:

template <typename T, index_t N, typename H>
struct Digest<Rect<T,N>, H> {
    H operator()(const Rect<T,N>& r) const {
        H nonce = geom::truncated_constant<H>(0x5fdb6f041e87e25d, 0x3ee84a498958343d);
        return geom::hash_many<H>(nonce, r.lo, r.hi);
    }
};

#ifdef GEOMC_USE_STREAMS

template <typename T, index_t N>
std::ostream& operator<<(std::ostream& os, const Rect<T,N>& r) {
    os << "[" << r.lo << " : " << r.hi << "]";
    return os;
}

#endif


} // end namespace geom


template <typename T, index_t N>
struct std::hash<geom::Rect<T,N>> {
    size_t operator()(const geom::Rect<T,N> &r) const {
        return geom::hash<geom::Rect<T,N>, size_t>(r);
    }
};
