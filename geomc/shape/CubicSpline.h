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

/**
 * @defgroup spline
 * @brief Splines and curves.
 */

/**
 * @addtogroup spline
 * @{
 */

/**
 * @brief Concept for a cubic spline.
 *
 * A cubic spline is a curve defined by four control points, and a basis matrix.
 * The control points have different meanings for each type of spline,
 * but all cubic splines are interconvertible.
 *
 * The basis matrix converts control points to cubic polynomial coefficients.
 * The resultant vector has the coefficients in order of increasing power.
 */
template <typename S, typename T, index_t N>
concept CubicSplineObject = requires (S s, const VecType<T,N> pts[4]) {
    S(pts);
    S(pts[0], pts[1], pts[2], pts[3]);
    { S::basis()         } -> std::convertible_to<const SimpleMatrix<T,4,4>&>;
    { S::inverse_basis() } -> std::convertible_to<const SimpleMatrix<T,4,4>&>;
    { s.control_points() } -> std::convertible_to<const VecType<T,N>*>;
    { s.bounds()         } -> std::same_as<Rect<T,N>>;
    { s.velocity(T(0))   } -> std::same_as<VecType<T,N>>;
    { s(T(0))            } -> std::same_as<VecType<T,N>>;
};


/**
 * @brief A cubic polynomial spline.
 *
 * PolynomialSplines store the coefficients of the polynomial directly.
 * They are more efficient to evaluate than other cubic splines, but are not
 * as convenient to work with.
 *
 * If evaluating a spline many times, it is more efficient to convert it to a
 * PolynomialSpline first.
 */
template <typename T, index_t N>
class PolynomialSpline : public Dimensional<T, N> {
public:
    /// Constant coefficient
    VecType<T,N> k0;
    /// Coefficient for s
    VecType<T,N> k1;
    /// Coefficient for s<sup>2</sup>.
    VecType<T,N> k2;
    /// Coefficient for s<sup>3</sup>.
    VecType<T,N> k3;
    
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
    constexpr PolynomialSpline(VecType<T,N> k0, VecType<T,N> k1, VecType<T,N> k2, VecType<T,N> k3):
        k0(k0), k1(k1), k2(k2), k3(k3) {}
    /// Construct a cubic polynomial spline from an array of coefficients.
    constexpr PolynomialSpline(const VecType<T,N> pts[4]):
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
    
    /// Get a const array of the four control points.
    constexpr const VecType<T,N>* control_points() const { return &k0; }
    /// Get an array of the four control points.
    constexpr       VecType<T,N>* control_points()       { return &k0; }
    
    /// The basis matrix for a cubic polynomial spline, which is the identity matrix.
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    /// The inverse basis matrix for a cubic polynomial spline, which is the identity matrix.
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
    /// Compute the bounding box of the spline.
    Rect<T,N> bounds() const {
        Rect<T,N> bounds = Rect<T,N>::from_corners(
            k0,               // s = 0 point
            k0 + k1 + k2 + k3 // s = 1 point
        );
        // look for inflection points between the endpoints.
        // where does each coordinate reach an extreme?
        for (index_t axis = 0; axis < N; ++axis) {
            T k = Vec<T,4> {
                coord(k0, axis),
                coord(k1, axis),
                coord(k2, axis),
                coord(k3, axis)
            };
            T results[2];
            bool ok = geom::quadratic_solve(results, k[3] * 3, k[2] * 2, k[1]);
            if (ok) {
                for (index_t i = 0; i < 2; ++i) {
                    T s = results[i];
                    if (s > 0 and s < 1) {
                        T s2 = s * s;
                        T v = k.dot({1, s, s2, s2 * s});
                        coord(bounds.lo, axis) = std::min(bounds.lo, v);
                        coord(bounds.hi, axis) = std::max(bounds.hi, v);
                    }
                }
            }
        }
        return bounds;
    }
    
    /// Evaluate the spline at a given parameter value.
    constexpr VecType<T,N> operator()(T s) const {
        T s2 = s * s;
        return k0 + k1 * s + k2 * s2 + k3 * s2 * s;
    }
    
    /// Compute the derivative (velocity) of the spline at a given parameter value.
    constexpr VecType<T,N> velocity(T s) const {
        return k1 + k2 * 2 * s + k3 * 3 * s * s;
    }
    
    /// Compute the second derivative (acceleration) of the spline at a given parameter value.
    constexpr VecType<T,N> acceleration(T s) const {
        return k2 * 2 + k3 * 6 * s;
    }
    
};

/**
 * @brief Base class for cubic splines.
 */
template <typename T, index_t N, typename Derived>
class CubicSpline : public Dimensional<T,N> {
    
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
          Derived& derived()       { return static_cast<      Derived&>(*this); }
    
    template <typename U, index_t M, typename S>
    S _make_spline(const SimpleMatrix<T,4,N>& mx) const {
        if constexpr (M > 1) {
            return S {
                mx.row(0),
                mx.row(1),
                mx.row(2),
                mx.row(3)
            };
        } else {
            return S {
                mx(0, 0),
                mx(1, 0),
                mx(2, 0),
                mx(3, 0)
            };
        }
    }
    
public:
    
    /// Convert to another type of cubic spline.
    template <CubicSplineObject<T,N> Spline>
    requires (not std::same_as<Spline, PolynomialSpline<T,N>>)
    constexpr operator Spline() const {
        using ptype = PointType<T,N>;
        WrapperMatrix<T,4,N> b = {ptype::iterator(derived().control_points()[0])};
        // (4 x 4) * (4 x N) = (4 x N)
        SimpleMatrix<T,4,N> coeffs = Spline::inverse_basis() * Derived::basis() * b;
        return _make_spline<T,N,Spline>(coeffs);
    }
    
    /// Convert this spline to its coefficient representation.
    constexpr operator PolynomialSpline<T,N>() const {
        using ptype = PointType<T,N>;
        WrapperMatrix<T,4,N> b = {ptype::iterator(derived().control_points()[0])};
        SimpleMatrix<T,4,N> coeffs = Derived::basis() * b;
        return _make_spline<T,N,PolynomialSpline<T,N>>(coeffs);
    }
    
    /// Evaluate the spline at a given parameter value.
    constexpr Vec<T,N> operator()(T s) const {
        PolynomialSpline<T,N> poly {*this};
        return poly(s);
    }
    
    /// Compute the derivative (velocity) of the spline at a given parameter value.
    constexpr Vec<T,N> velocity(T s) const {
        PolynomialSpline<T,N> poly {*this};
        return poly.velocity(s);
    }
    
    /// Compute the second derivative (acceleration) of the spline at a given parameter value.
    constexpr Vec<T,N> acceleration(T s) const {
        PolynomialSpline<T,N> poly {*this};
        return poly.acceleration(s);
    }
    
    /// Compute the bounding box of the spline.
    constexpr Rect<T,N> bounds() const {
        PolynomialSpline<T,N> poly {*this};
        return poly.bounds();
    }
};


/**
 * @brief A cubic spline with two knots and two tangent points.
 */
template <typename T, index_t N>
class BezierSpline : public CubicSpline<T,N,BezierSpline<T,N>> {
public:
    /// Position at s = 0
    VecType<T,N> p0;
    /// Tangent point for s = 0
    VecType<T,N> t0;
    /// Tangent point for s = 1
    VecType<T,N> t1;
    /// Position at s = 1
    VecType<T,N> p1;
    
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
    constexpr BezierSpline(VecType<T,N> p0, VecType<T,N> t0, VecType<T,N> t1, VecType<T,N> p1):
        p0(p0), t0(t0), t1(t1), p1(p1) {}
    /// Construct the spline from an array of control points.
    constexpr BezierSpline(const VecType<T,N> pts[4]):
        p0(pts[0]), t0(pts[1]), t1(pts[2]), p1(pts[3]) {}
    
    /// Get the control points.
    constexpr const VecType<T,N>* control_points() const { return &p0; }
    /// Get the control points.
    constexpr       VecType<T,N>* control_points()       { return &p0; }
    
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
 * Transform a Bezier spline.
 * @related BezierSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
BezierSpline<T,N> operator*(
    const Xf& xf,
    const BezierSpline<T,N>& spline)
{
    return {
        xf * spline.p0,
        xf * spline.t0,
        xf * spline.t1,
        xf * spline.p1
    };
}

/**
 * Inverse transform a Bezier spline.
 * @related BezierSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
BezierSpline<T,N> operator/(
    const BezierSpline<T,N>& spline,
    const Xf& xf)
{
    return {
        spline.p0 / xf,
        spline.t0 / xf,
        spline.t1 / xf,
        spline.p1 / xf
    };
}


/**
 * @brief A cubic spline with four guide points and continuous curvature.
 *
 * The curve does not necessarily pass through the guide points.
 * Concatenated B-splines which share three consecutive guide points are
 * continuous in position, velocity, and curvature.
 */
template <typename T, index_t N>
class BSpline : public CubicSpline<T,N,BSpline<T,N>> {
public:
    /// Guide point for s = -1
    VecType<T,N> s0;
    /// Guide point for s = 0
    VecType<T,N> s1;
    /// Guide point for s = 1
    VecType<T,N> s2;
    /// Guide point for s = 2
    VecType<T,N> s3;
    
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
    constexpr BSpline(VecType<T,N> s0, VecType<T,N> s1, VecType<T,N> s2, VecType<T,N> s3):
        s0(s0), s1(s1), s2(s2), s3(s3) {}
    /// Construct the spline from an array of control points.
    constexpr BSpline(const VecType<T,N> pts[4]):
        s0(pts[0]), s1(pts[1]), s2(pts[2]), s3(pts[3]) {}
    
    constexpr const VecType<T,N>* control_points() const { return &s0; }
    constexpr       VecType<T,N>* control_points()       { return &s0; }
    
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
};

/**
 * Transform a B-spline.
 * @related BSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
BSpline<T,N> operator*(
    const Xf& xf,
    const BSpline<T,N>& spline)
{
    return {
        xf * spline.s0,
        xf * spline.s1,
        xf * spline.s2,
        xf * spline.s3
    };
}

/**
 * Inverse transform a B-spline.
 * @related BSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
BSpline<T,N> operator/(
    const BSpline<T,N>& spline,
    const Xf& xf)
{
    return {
        spline.s0 / xf,
        spline.s1 / xf,
        spline.s2 / xf,
        spline.s3 / xf
    };
}


/**
 * @brief A cubic spline which passes smoothly through four knots.
 *
 * Concatenated Catmull-Rom splines are continuous in position and velocity at the
 * control points, but not in curvature.
 */
template <typename T, index_t N>
class CatromSpline : public CubicSpline<T,N,CatromSpline<T,N>> {
public:
    /// Position at s = -1
    VecType<T,N> p0;
    /// Position at s = 0
    VecType<T,N> p1;
    /// Position at s = 1
    VecType<T,N> p2;
    /// Position at s = 2
    VecType<T,N> p3;
    
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
    constexpr CatromSpline(VecType<T,N> p0, VecType<T,N> p1, VecType<T,N> p2, VecType<T,N> p3):
        p0(p0), p1(p1), p2(p2), p3(p3) {}
    /// Construct the spline from an array of control points.
    constexpr CatromSpline(const VecType<T,N> pts[4]):
        p0(pts[0]), p1(pts[1]), p2(pts[2]), p3(pts[3]) {}
    
    constexpr const VecType<T,N>* control_points() const { return &p0; }
    constexpr       VecType<T,N>* control_points()       { return &p0; }
    
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
    /// Compute a transformation of the spline.
    template <Transform<T,N> Xf>
    CatromSpline operator*(Xf xf) const {
        return BezierSpline {
            xf * p0,
            xf * p1,
            xf * p2,
            xf * p3
        };
    }
    
};

/**
 * Transform a Catmull-Rom spline.
 * @related CatromSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
CatromSpline<T,N> operator*(
    const Xf& xf,
    const CatromSpline<T,N>& spline)
{
    return {
        xf * spline.p0,
        xf * spline.p1,
        xf * spline.p2,
        xf * spline.p3
    };
}

/**
 * Inverse transform a Catmull-Rom spline.
 * @related CatromSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
CatromSpline<T,N> operator/(
    const CatromSpline<T,N>& spline,
    const Xf& xf)
{
    return {
        spline.p0 / xf,
        spline.p1 / xf,
        spline.p2 / xf,
        spline.p3 / xf
    };
}


/**
 * @brief A cubic spline defined by two points and two velocities.
 */
template <typename T, index_t N>
class HermiteSpline : public CubicSpline<T,N,HermiteSpline<T,N>> {
public:
    /// Position at s = 0
    VecType<T,N> p0;
    /// Velocity at s = 0
    VecType<T,N> v0;
    /// Position at s = 1
    VecType<T,N> p1;
    /// Velocity at s = 1
    VecType<T,N> v1;
    
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
    constexpr HermiteSpline(VecType<T,N> p0, VecType<T,N> v0, VecType<T,N> p1, VecType<T,N> v1):
        p0(p0), v0(v0), p1(p1), v1(v1) {}
    /// Construct the spline from an array of control points.
    constexpr HermiteSpline(const VecType<T,N> pts[4]):
        p0(pts[0]), v0(pts[1]), p1(pts[2]), v1(pts[3]) {}
    
    constexpr const VecType<T,N>* control_points() const { return &p0; }
    constexpr       VecType<T,N>* control_points()       { return &p0; }
    
    static constexpr SimpleMatrix<T,4,4> basis() {
        return ControlToCoefficients;
    }
    
    static constexpr SimpleMatrix<T,4,4> inverse_basis() {
        return CoefficientsToControl;
    }
    
};

/**
 * Transform a Hermite spline.
 * @related HermiteSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
HermiteSpline<T,N> operator*(
    const Xf& xf,
    const HermiteSpline<T,N>& spline)
{
    return {
        xf * spline.p0,
        xf.apply_direction(spline.v0),
        xf * spline.p1,
        xf.apply_direction(spline.v1)
    };
}

/**
 * Inverse transform a Hermite spline.
 * @related HermiteSpline
 */
template <typename T, index_t N, Transform<T,N> Xf>
HermiteSpline<T,N> operator/(
    const HermiteSpline<T,N>& spline,
    const Xf& xf)
{
    return {
        spline.p0 / xf,
        xf.apply_inverse_direction(spline.v0),
        spline.p1 / xf,
        xf.apply_inverse_direction(spline.v1)
    };
}


/// @} // addtogroup spline
/// @} // addtogroup shape

} // namespace geom
