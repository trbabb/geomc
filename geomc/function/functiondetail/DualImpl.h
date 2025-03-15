#pragma once

#include <geomc/geomc_defs.h>


namespace geom {

/**
 * @brief Policy for handling discontinuous derivatives
 *
 * @ingroup function
 */
enum class DiscontinuityPolicy {
    /// Return `NaN` at discontinuities
    NaN,
    /// At discontinuities, return the average of the two boundary values.
    Average,
    /// At discontinuities, return the value when approaching from the left
    Left,
    /// At discontinuities, return the value when approaching from the right
    Right,
    /**
     * @brief At discontinuities, return +/- infinity, according to whether the
     * discontinuity is an increasing or decreasing jump.
     */
    Inf
};

namespace detail {

template <typename T, DiscontinuityPolicy P>
struct discontinuity_impl {};

template <typename T>
struct discontinuity_impl<T, DiscontinuityPolicy::NaN> {
    static inline T resolve(T a, T b) {
        return std::numeric_limits<T>::quiet_NaN();
    }
};

template <typename T>
struct discontinuity_impl<T, DiscontinuityPolicy::Left> {
    static inline T resolve(T a, T b) {
        return a;
    }
};

template <typename T>
struct discontinuity_impl<T, DiscontinuityPolicy::Right> {
    static inline T resolve(T a, T b) {
        return b;
    }
};

template <typename T>
struct discontinuity_impl<T, DiscontinuityPolicy::Average> {
    static inline T resolve(T a, T b) {
        return (a + b) / (T)2;
    }
};

template <typename T>
struct discontinuity_impl<T, DiscontinuityPolicy::Inf> {
    static inline T resolve(T a, T b) {
        const T inf = std::numeric_limits<T>::infinity();
        return (b == a) ? a : (b > a ? inf : -inf);
    }
};

} // namespace detail


/**
 * @brief Class implementing the dual numbers, whose arithmetic operations perform a 
 * simultaneous calculation of the first derivative.
 * 
 * Overview
 * ========
 * 
 * Dual numbers implicitly compute and store the first derivative of any arbitrary 
 * function by extending the reals to include a new nonzero element &epsilon; 
 * whose square is zero. Thus duals have the form (a + b&epsilon;); where `a` 
 * may be thought of as the function value, and `b` its first derivative at `a`. 
 * This technique is commonly known as Automatic Differentiation, or AD.
 * 
 * Fundamental operations on duals (including addition, multiplication, division, 
 * and any number of other primitive mathematical functions) implicitly compute and 
 * keep track of the derivative by performing the chain and product rules in-place.
 * Efficiency and accuracy is very good (easily better than either symbolic or 
 * numerical differentiation), adding only a constant factor to arithmetic
 * operations.
 * 
 * In general, Duals should behave exactly as bare numbers in code.
 *
 * It is safe to convert a constant value to a Dual (and this will be done implicitly
 * by setting the epsilon component to zero). It is unsafe to convert a Dual number
 * to a bare value, as this forgets the derivative information. Therefore if the
 * Dual is to be used as a constant, it must be cast explicitly, or the primal
 * component must be extracted.
 * 
 * Use
 * ===
 * 
 * The value of an arbitrary function may be differentiated with respect
 * to its input if and only if those inputs, and all the calculations that depend on
 * those inputs, are performed entirely on Dual numbers. In other words, a function
 * is differentiable if it uses Duals "all the way to the bottom" of its calculation.
 * If ever the primal is extracted and used in some sub-calculation apart from its
 * epsilon component, then the resultant derivative information will be invalid.
 * 
 * The epsilon components of a function's input or output together represent the
 * components of a directional derivative.
 * 
 * Equality
 * ========
 *
 * Duals are compared only using their primal value, so that they behave just as
 * bare numbers do. This means that two `Dual`s with different epsilons may compare
 * equal! If you need exact binary equality, then explicitly check that both `x`s
 * and `dx`s are equal. 
 *
 * @tparam T The numerical type of the number.
 * @tparam DP The `DiscontinuityPolicy` for handling discontinuous derivatives.
 * 
 * @ingroup function
 */
template <typename T, DiscontinuityPolicy Dp=DiscontinuityPolicy::Right>
struct Dual {
    typedef detail::discontinuity_impl<T, Dp> discontinuity_policy;
    
    /// Real (primal) component
    T x;
    /// Dual (epsilon) component
    T dx;
    
    /*******************************
     * Constructors                *
     *******************************/
   
    /// Construct a Dual default value and 0 derivative.
    constexpr Dual():x(0), dx(0) {}
    
    /// Construct a Dual with value `x` and 0 derivative.
    constexpr Dual(T x):x(x), dx(0) {}
    
    /// Construct a Dual with value `x` and derivative `dx`.
    constexpr Dual(T x, T dx):x(x), dx(dx) {}
    
    /// Construct a Dual from one with a different discontinuity policy.
    template <DiscontinuityPolicy DP1>
    constexpr Dual(const Dual<T,DP1>& d):x(d.dx), dx(d.dx) {}
    
    /// Copy-construct a Dual
    constexpr Dual(const Dual& d) = default;
    
    
    /*******************************
     * Operators                   *
     *******************************/
    
    // assign
    
    constexpr Dual<T,Dp>& operator=(const Dual<T,Dp>& other) = default;
    constexpr Dual<T,Dp>& operator=(T other) {
        x  = other;
        dx = 0;
        return *this;
    }
    
    // comparison
    
    /// Equality of primal component
    inline bool operator==(const Dual<T,Dp> &d) { return x == d.x; }
    /// Inequality of primal component
    inline bool operator!=(const Dual<T,Dp> &d) { return x != d.x; }
    /// Greater-or-equal-to comparison of primal component
    inline bool operator>=(const Dual<T,Dp> &d) { return x >= d.x; }
    /// Less-or-equal-to comparison of primal component
    inline bool operator<=(const Dual<T,Dp> &d) { return x <= d.x; }
    /// Greater-than comparison of primal component
    inline bool operator> (const Dual<T,Dp> &d) { return x  > d.x; }
    /// Less-than comparison of primal component
    inline bool operator< (const Dual<T,Dp> &d) { return x  < d.x; }
    
    /// Post-increment
    Dual<T,Dp> operator++(int) {
        Dual<T,Dp> d = *this;
        x++;
        return d;
    }
    
    /// Post-decrement
    Dual<T,Dp> operator--(int) {
        Dual<T,Dp> d = *this;
        x--;
        return d;
    }
    
    /// Pre-increment
    Dual<T,Dp>& operator++() {
        x++;
        return *this;
    }
    
    /// Pre-decrement
    Dual<T,Dp>& operator--() {
        x--;
        return *this;
    }
    
    // conversion
    
    template <typename U>
    requires std::convertible_to<T,U>
    inline operator Dual<U>() const {
        return {
            static_cast<U>(x), 
            static_cast<U>(dx)
        };
    }
    
    template <typename U, DiscontinuityPolicy P2=Dp>
    requires std::convertible_to<T,U>
    explicit inline operator U() const {
        return static_cast<U>(x);
    }

}; // struct Dual


#if GEOMC_USE_STREAMS

// xxx todo: this isn't found. why?
template <typename T, DiscontinuityPolicy Dp>
std::ostream& operator<<(std::ostream& out, const Dual<T,Dp>& d) {
    out << d.x << " + " << d.dx << "dx";
    return out;
}

#endif

} // namespace geom
