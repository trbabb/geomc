#pragma once

#include <geomc/function/Utils.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Matrix.h>
#include <geomc/shape/ShapeTypes.h>

// todo: make sure this works with N = 1

namespace geom {

/**
 * @addtogroup shape
 * @{
 */

template <typename T, index_t N>
constexpr Vec<T,N> eval_cubic_spline(
    const SimpleMatrix<T,4,4>& basis,
    const VecType<T,N> control_points[4],
    T s)
{
    T s2 = s * s;
    Vec<T,4> s_pow = {1, s, s2, s2 * s};
    WrapperMatrix<T,4,N> b = {control_points};
    // (4 x 4) * (4 x N) = (4 x N)
    SimpleMatrix<T,4,N> coeffs = basis * b;
    // (1 x 4) * (4 x N) = (1 x N)
    return s_pow * coeffs;
}

template <typename T, index_t N>
constexpr Vec<T,N> eval_cubic_spline_derivative(
    const SimpleMatrix<T,4,4>& basis,
    const VecType<T,N> control_points[4],
    T s)
{
    Vec<T,4> s_pow = {0, 1, 2 * s, 3 * s * s};
    WrapperMatrix<T,4,N> b = {control_points};
    // (4 x 4) * (4 x N) = (4 x N)
    SimpleMatrix<T,4,N> coeffs = basis * b;
    // (1 x 4) * (4 x N) = (1 x N)
    return s_pow * coeffs;
}

template <typename T, index_t N>
constexpr Rect<T,N> cubic_spline_bounds(
    const SimpleMatrix<T,4,4>& basis,
    const VecType<T,N> control_points[4])
{
    Rect<T,N> bounds;
    WrapperMatrix<T,4,N> b = {control_points};
    // (4 x 4) * (4 x N) = (4 x N)
    SimpleMatrix<T,4,N> coeffs = basis * b;
    bounds |= coeffs.row(0);        // s = 0 point
    bounds |= Vec<T,4>(1) * coeffs; // s = 1 point
    // look for inflection points between the endpoints.
    // where does each coordinate reach an extreme?
    for (index_t axis = 0; axis < N; ++axis) {
        T a = coeffs(3, axis) * 3;
        T b = coeffs(2, axis) * 2;
        T c = coeffs(1, axis);
        T results[2];
        bool ok = geom::quadratic_solve(results, a, b, c);
        if (ok) {
            for (index_t i = 0; i < 2; ++i) {
                T s = results[i];
                if (s > 0 and s < 1) {
                    T s2 = s * s;
                    T v = Vec<T,N>{1, s, s2, s2 * s} * coeffs.col(axis);
                    bounds.lo = std::min(bounds.lo, v);
                    bounds.hi = std::max(bounds.hi, v);
                }
            }
        }
    }
    return bounds;
}

/**
 * @brief Concept for a cubic spline.
 *
 * A cubic spline is a curve defined by four control points, and a basis matrix.
 * The control points have different meanings for each type of spline,
 * but all cubic splines are interconvertible.
 */
template <typename S, typename T, index_t N>
concept CubicSplineObject = requires (S s, const VecType<T,N> pts[4]) {
    S(pts);
    S(pts[0], pts[1], pts[2], pts[3]);
    { s.control_points() } -> std::convertible_to<const VecType<T,N>*>;
    { s.basis()          } -> std::convertible_to<const SimpleMatrix<T,4,4>&>;
    { s.inverse_basis()  } -> std::convertible_to<const SimpleMatrix<T,4,4>&>;
    { s.bounds()         } -> std::same_as<Rect<T,N>>;
    { s.d_ds(T(0))       } -> std::same_as<VecType<T,N>>;
    { s(T(0))            } -> std::same_as<VecType<T,N>>;
};


/**
 * @brief Base class for cubic splines.
 */
template <typename T, index_t N, typename Derived>
class CubicSpline : public Dimensional<T,N> {
    
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
          Derived& derived()       { return static_cast<      Derived&>(*this); }
    
public:
    
    /// Convert to another type of cubic spline.
    template <CubicSplineObject<T,N> Spline>
    requires (not std::same_as<Spline, PolynomialSpline<T,N>>)
    constexpr operator Spline() const {
        WrapperMatrix<T,4,N> b = {derived().control_points()};
        // (4 x 4) * (4 x N) = (4 x N)
        SimpleMatrix<T,4,N> coeffs = Spline::inverse_basis() * Derived::basis() * b;
        return Spline {
            coeffs.row(0),
            coeffs.row(1),
            coeffs.row(2),
            coeffs.row(3)
        };
    }
    
    /// Convert this spline to its coefficient representation.
    constexpr operator PolynomialSpline<T,N>() const {
        WrapperMatrix<T,4,N> b = {derived().control_points()};
        SimpleMatrix<T,4,N> coeffs = Derived::basis() * b;
        return PolynomialSpline<T,N> {
            coeffs.row(0),
            coeffs.row(1),
            coeffs.row(2),
            coeffs.row(3)
        };
    }
    
    /// Evaluate the spline at a given parameter value.
    constexpr Vec<T,N> operator()(T s) const {
        eval_cubic_spline(Derived::basis(), derived().control_points(), s);
    }
    
    /// Compute the derivative (velocity) of the spline at a given parameter value.
    constexpr Vec<T,N> d_ds(T s) const {
        eval_cubic_spline_derivative(Derived::basis(), derived().control_points(), s);
    }
    
    /// Compute the bounding box of the spline.
    constexpr Rect<T,N> bounds() const {
        cubic_spline_bounds(Derived::basis(), derived().control_points());
    }
};


/**
 * @brief A cubic Bezier curve.
 *
 * A Bezier curve is defined by four control points: the start and end points,
 * and the tangent points at the start and end.
 */
template <typename T, index_t N>
class BezierSpline : public CubicSpline<T,N,BezierSpline<T,N>> {
public:
    /// Position at s = 0
    Vec<T,N> p0;
    /// Tangent point for s = 0
    Vec<T,N> t0;
    /// Tangent point for s = 1
    Vec<T,N> t1;
    /// Position at s = 1
    Vec<T,N> p1;
    
    static constexpr SimpleMatrix<T,4,4> ControlToCoefficients = {
        { 1,  0,  0, 0},
        {-3,  3,  0, 0},
        { 3, -6,  3, 0},
        {-1,  3, -3, 1},
    };
    
    static constexpr SimpleMatrix<T,4,4> CoefficientsToControl = {
        {3, 0, 0, 0},
        {3, 1, 0, 0},
        {3, 2, 1, 0},
        {3, 3, 3, 3},
    };
    
    constexpr BezierSpline() {
        p1[0] = 1;
        t0[0] = 0.25;
        t1[0] = 0.75;
    }
    
    /// Construct the spline from four control points.
    constexpr BezierSpline(Vec<T,N> p0, Vec<T,N> t0, Vec<T,N> t1, Vec<T,N> p1):
        p0(p0), t0(t0), t1(t1), p1(p1) {}
    /// Construct the spline from an array of control points.
    constexpr BezierSpline(const Vec<T,N> pts[4]):
        p0(pts[0]), t0(pts[1]), t1(pts[2]), p1(pts[3]) {}
    
    /// Get the control points.
    constexpr const Vec<T,N>* control_points() const { return &p0; }
    /// Get the control points.
    constexpr       Vec<T,N>* control_points()       { return &p0; }
    
    /// Get the basis matrix, which converts control points to coefficients.
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    /// Get the inverse basis matrix, which converts coefficients to control points.
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
};


/**
 * @brief A cubic B-spline.
 *
 * B-spline curves have continuous curvature at the control points.
 *
 * The control points guide the shape of the curve, but the curve does not pass through them
 * in general.
 */
template <typename T, index_t N>
class BSpline : public CubicSpline<T,N,BSpline<T,N>> {
public:
    /// Guide point for s = -1
    Vec<T,N> s0;
    /// Guide point for s = 0
    Vec<T,N> s1;
    /// Guide point for s = 1
    Vec<T,N> s2;
    /// Guide point for s = 2
    Vec<T,N> s3;
    
    static constexpr SimpleMatrix<T,4,4> ControlToCoefficients = {
        { 1,  4,  1, 0},
        {-3,  0,  3, 0},
        { 3, -6,  3, 0},
        {-1,  3, -3, 1},
    };
    
    static constexpr SimpleMatrix<T,4,4> CoefficientsToControl = {
        {3, -3,  2,  0},
        {3,  0, -1,  0},
        {3,  3,  2,  0},
        {3,  6, 11, 18},
    };
    
    constexpr BSpline() {
        s0[0] = -1;
        s1[0] =  0;
        s2[0] =  1;
        s3[0] =  2;
    }
    
    /// Construct the spline from four control points.
    constexpr BSpline(Vec<T,N> s0, Vec<T,N> s1, Vec<T,N> s2, Vec<T,N> s3):
        s0(s0), s1(s1), s2(s2), s3(s3) {}
    /// Construct the spline from an array of control points.
    constexpr BSpline(const Vec<T,N> pts[4]):
        s0(pts[0]), s1(pts[1]), s2(pts[2]), s3(pts[3]) {}
    
    constexpr const Vec<T,N>* control_points() const { return &s0; }
    constexpr       Vec<T,N>* control_points()       { return &s0; }
    
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
};


/**
 * @brief A cubic Catmull-Rom spline.
 *
 * Catmull-Rom splines are defined by four control points, which the curve passes through.
 * The curve is continuous in position and velocity at the control points, but not in
 * curvature.
 */
template <typename T, index_t N>
class CatromSpline : public CubicSpline<T,N,CatromSpline<T,N>> {
public:
    /// Position at s = -1
    Vec<T,N> p0;
    /// Position at s = 0
    Vec<T,N> p1;
    /// Position at s = 1
    Vec<T,N> p2;
    /// Position at s = 2
    Vec<T,N> p3;
    
    static constexpr SimpleMatrix<T,4,4> ControlToCoefficients = {
        { 0,  2,  0,  0},
        {-1,  0,  1,  0},
        { 2, -5,  4, -1},
        {-1,  3, -3,  1},
    };
    
    static constexpr SimpleMatrix<T,4,4> CoefficientsToControl = {
        {-1, -3, -1, -1},
        { 0, -1, -2, -1},
        {-3,  4, -5,  0},
        { 0, -4,  2, -2},
    };
    
    constexpr CatromSpline() {
        p0[0] = -1;
        p1[0] =  0;
        p2[0] =  1;
        p3[0] =  2;
    }
    
    /// Construct the spline from four control points.
    constexpr CatromSpline(Vec<T,N> p0, Vec<T,N> p1, Vec<T,N> p2, Vec<T,N> p3):
        p0(p0), p1(p1), p2(p2), p3(p3) {}
    /// Construct the spline from an array of control points.
    constexpr CatromSpline(const Vec<T,N> pts[4]):
        p0(pts[0]), p1(pts[1]), p2(pts[2]), p3(pts[3]) {}
    
    constexpr const Vec<T,N>* control_points() const { return &p0; }
    constexpr       Vec<T,N>* control_points()       { return &p0; }
    
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
};


/**
 * @brief A cubic Hermite spline.
 *
 * Hermite splines are defined by four control points: the start and end points,
 * and the velocities at the start and end.
 */
template <typename T, index_t N>
class HermiteSpline : public CubicSpline<T,N,HermiteSpline<T,N>> {
public:
    /// Position at s = 0
    Vec<T,N> p0;
    /// Velocity at s = 0
    Vec<T,N> v0;
    /// Position at s = 1
    Vec<T,N> p1;
    /// Velocity at s = 1
    Vec<T,N> v1;
    
    static constexpr SimpleMatrix<T,4,4> ControlToCoefficients = {
        { 1,  0,  0,  0},
        { 0,  1,  0,  0},
        {-3, -2,  3, -1},
        { 2,  1, -2,  1},
    };
    
    static constexpr SimpleMatrix<T,4,4> CoefficientsToControl = {
        {-2, -1, -1, -1},
        {-1, -2, -1, -1},
        { 2,  1, -4,  0},
        {-3, -2,  1, -2},
    };
    
    constexpr HermiteSpline() {
        v0[0] = 1;
        p1[0] = 1;
        v1[0] = 1;
    }
    
    /// Construct the spline from four control points.
    constexpr HermiteSpline(Vec<T,N> p0, Vec<T,N> v0, Vec<T,N> p1, Vec<T,N> v1):
        p0(p0), v0(v0), p1(p1), v1(v1) {}
    /// Construct the spline from an array of control points.
    constexpr HermiteSpline(const Vec<T,N> pts[4]):
        p0(pts[0]), v0(pts[1]), p1(pts[2]), v1(pts[3]) {}
    
    constexpr const Vec<T,N>* control_points() const { return &p0; }
    constexpr       Vec<T,N>* control_points()       { return &p0; }
    
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
};


/**
 * @brief A cubic polynomial spline.
 *
 * PolynomialSplines store the coefficients of the polynomial directly.
 * They are more efficient to evaluate than other cubic splines, but are not
 * as convenient to work with.
 */
template <typename T, index_t N>
class PolynomialSpline : public CubicSpline<T,N,PolynomialSpline<T,N>> {
public:
    /// Constant coefficient
    Vec<T,N> k0;
    /// Coefficient for s
    Vec<T,N> k1;
    /// Coefficient for s^2
    Vec<T,N> k2;
    /// Coefficient for s^3
    Vec<T,N> k3;
    
    static constexpr SimpleMatrix<T,4,4> ControlToCoefficients = {
        { 1,  0,  0,  0},
        { 0,  1,  0,  0},
        { 0,  0,  1,  0},
        { 0,  0,  0,  1},
    };
    
    static constexpr SimpleMatrix<T,4,4> CoefficientsToControl = {
        { 1,  0,  0,  0},
        { 0,  1,  0,  0},
        { 0,  0,  1,  0},
        { 0,  0,  0,  1},
    };
    
    constexpr PolynomialSpline() {
        k1[0] = 1;
    }
    
    /// Construct a cubic polynomial spline from coefficients.
    constexpr PolynomialSpline(Vec<T,N> k0, Vec<T,N> k1, Vec<T,N> k2, Vec<T,N> k3):
        k0(k0), k1(k1), k2(k2), k3(k3) {}
    /// Construct a cubic polynomial spline from an array of coefficients.
    constexpr PolynomialSpline(const Vec<T,N> pts[4]):
        k0(pts[0]), k1(pts[1]), k2(pts[2]), k3(pts[3]) {}
    
    /// Convert to another type of cubic spline.
    template <CubicSplineObject<T,N> Spline>
    constexpr operator Spline() const {
        WrapperMatrix<T,4,N> b = {control_points()};
        // (4 x 4) * (4 x N) = (4 x N)
        SimpleMatrix<T,4,N> coeffs = Spline::inverse_basis() * b;
        return Spline {
            coeffs.row(0),
            coeffs.row(1),
            coeffs.row(2),
            coeffs.row(3)
        };
    }
    
    constexpr const Vec<T,N>* control_points() const { return &k0; }
    constexpr       Vec<T,N>* control_points()       { return &k0; }
    
    /// The basis matrix for a cubic polynomial spline, which is the identity matrix.
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    /// The inverse basis matrix for a cubic polynomial spline, which is the identity matrix.
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
    /// Evaluate the spline at a given parameter value.
    constexpr Vec<T,N> operator()(T s) const {
        T s2 = s * s;
        return k0 + k1 * s + k2 * s2 + k3 * s2 * s;
    }
    
    /// Compute the derivative (velocity) of the spline at a given parameter value.
    constexpr Vec<T,N> d_ds(T s) const {
        return k1 + k2 * 2 * s + k3 * 3 * s * s;
    }
    
};

/// @} // addtogroup shape

} // namespace geom
