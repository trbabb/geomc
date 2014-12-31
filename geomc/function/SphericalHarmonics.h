/* 
 * File:   SphericalHarmonics.h
 * Author: tbabb
 *
 * Created on December 26, 2014, 1:51 AM
 */

#ifndef SPHERICALHARMONICS_H
#define	SPHERICALHARMONICS_H

#include <cmath>
#include <math.h> // xxx debug
#include <geomc/Storage.h>
#include <geomc/function/Basis.h>

namespace geom {
    
/*
 * the final problem is how to organize this.
 * - how to deal with reconstruction
 *   it's a process of giving n samples and dotting them with the basis
 *   functions. we could do this easily if the n samples are in an array,
 *   but what if we want to do some online algorithm with a huge number of samples?
 *   or dynamic importance sampling/gradient descent/etc?
 * - how to handle euclidean samples? I just don't like the idea of trig on 
 *   the results of atan().
 *   - there exist expressions for this that involve combinatorics
 *     they are not formulated as recurrence relations, though you could
 *     probably do that work.
 *   - euclidean basis functions might allow N-dimensional SH.
 * 
 */

/* questions
 * 
 * - does the boost impl use condon-shortley phase?
 * - does the boost impl accept negative m?
 * - what's the deal with the sqrt(2) factor? Why does it only appear in some
 *   definitions of the SH basis functions?
 * - what is the proper handling of -m bands? Why do some definitions which
 *   state P(l,-m,x) in terms of P(l,m,x) include an extra partial factorial?
 *   why do these not appear in SH basis definitions that pass |m| to P(...)?
 *  
 */  

// todo: - verify basis functions
//         (reconstructions do not seem accurate)
//       - there is a 1-x^2 which can be computed from the angles, which are known.
//       - zonal harmonics/spot function
//       - (rotated) zonal harmonics -> SH
//       - sh rotation
//       - complex numbers
//       - convolution
//       - consider a `highest_nonzero` to cheapen
//         band-limiting ops?

/**
 * @addtogroup function
 * @{
 */

/**
 * Evaluate a spherical harmonic coefficient.
 * 
 * @code#include <geomc/function/SphericalHarmonics.h>@endcode
 * @param l Band index.
 * @param m Sub-band.
 * @param cos_alt Cosine of angle from polar axis.
 * @param azi Azimuthal (longitude) angle.
 */
template <typename T>
T spherical_harmonic_coeff(index_t l, index_t m, T cos_alt, T azi) {
    index_t m_a = std::abs(m);
    static const T sr2 = std::sqrt((T)2);
    static const T pi  = boost::math::constants::pi<T>();
    
    T P_lm = legendre(l, m_a, cos_alt);
    
    // k contains a sub-expression commonly written: (l-m)!/(l+m)!
    // which is vulnerable to precision/overflow issues, as very large
    // intermediate numbers must be stored. here, we compute the reduced
    // fraction directly:
    T z = 1;
    if (l != 0) for (index_t i = l - m_a + 1; i <= l + m_a; i++) { z *= i; }
    
    T k = std::sqrt((1/z) * ((2 * l + 1) / (4 * pi)));
    
    if (m > 0) {
        return sr2 * k * std::cos( m * azi) * P_lm;
    } else if (m < 0) {
        return sr2 * k * std::sin(-m * azi) * P_lm;
    } else {
        return k * legendre(l, 0, cos_alt);
    }
}

/**
 * @brief Class and methods for representing band-limited functions on the surface of a 
 * 3D sphere. 
 * 
 * @tparam T Type of evaluated function.
 * @tparam Bands Number of bands to represent, or 0 if the band count is dynamic. 
 * 
 * Memory requirements are O(n<sup>2</sup>) on the number of bands.
 * 
 * Complex-valued SH functions are not yet natively supported.
 */
template <typename T, index_t Bands>
class SphericalHarmonics {
    protected:
    typename Dimension<Bands>::storage_t n_bands;
    public:
    /// Coefficient vector.
    Storage<T,Bands*Bands> coeffs;
    
    /**
     * @brief Construct a SphericalHarmonics basis for the direction given by `d`.
     * @param d Direction for basis.
     * @param n Number of bands to use, if dynamic.
     */
    static SphericalHarmonics<T,Bands> basis(Vec<T,3> d, index_t n=Bands) {
        SphericalHarmonics<T,Bands> out(n);
        d = d.unit();
        out.project(d.z, std::atan2(d.y, d.x), 1);
        return out;
    }
    
    /**
     * Construct a new SphericalHarmonics. 
     * 
     * If dynamically sized, only a single DC band will be allocated.
     */
    SphericalHarmonics():coeffs(       std::max(Bands*Bands,(index_t)1)) {
        Dimension<Bands>::set(n_bands, std::max(Bands,      (index_t)1)); 
        std::fill(coeffs.get(), coeffs.get() + size(), 0);
    }
    
    /**
     * Construct a new SphericalHarmonics. (Dynamic size).
     * 
     * @param bands Number of bands to represent. Ignored unless the template 
     * parameter `Bands` is 0 (dynamic).
     */
    explicit SphericalHarmonics(index_t bands):coeffs(bands * bands) {
        Dimension<Bands>::set(n_bands, bands);
        std::fill(coeffs.get(), coeffs.get() + size(), 0);
    }
    
    /**
     * Get the stored SH coefficient for band `l` and sub-band `m`. 
     * @param l Integer in `[0, bands())`.
     * @param m Integer in `[-l, l]`.
     */
    T coeff(index_t l, index_t m) const {
        return coeffs[l * l + l + m];
    }
    
    /**
     * Get the stored SH coefficient for band `l` and sub-band `m`. 
     * @param l Integer in `[0, bands())`.
     * @param m Integer in `[-l, l]`.
     */
    T& coeff(index_t l, index_t m) {
        return coeffs[l * l + l + m];
    }
    
    /**
     * Get a pointer to the array of coefficients.
     */
    inline T* get() {
        return coeffs.get();
    }
    
    /**
     * Get a pointer to the array of coefficients.
     */
    inline const T* get() const {
        return coeffs.get();
    }
    
    /**
     * @return Number of coefficients stored.
     */
    inline index_t size() const {
        const index_t n = bands();
        return n * n;
    }
    
    /**
     * @return Number of bands stored.
     */
    inline index_t bands() const {
        return Dimension<Bands>::get(n_bands);
    }
    
    /**
     * Compute the inner product of this SH function with `other`.
     * @return Sum of products of respective coefficients between `this` and `other`.
     */
    T dot(const SphericalHarmonics<T,Bands> &other) const {
        const index_t n = std::min(size(), other.size());
        T sum = 0;
        for (index_t i = 0; i < n; i++) {
            sum += coeffs[i] * other.coeffs[i];
        }
        return sum;
    }
    
    
    protected:
        
    // protected implementation takes cos(alt), which is just z
    // for a unit vector. this formulation may better preserve precision when
    // z is available directly.
    T _eval(T cos_alt, T azi) const {
        static const T sr2 = std::sqrt((T)2);
        static const T pi  = boost::math::constants::pi<T>();
        const index_t l = bands();
        const T x = cos_alt;
        T sum = 0;
        for (index_t m = 0; m < l; m++) {
            T c_ma = std::cos(m * azi);
            T s_ma = std::sin(m * azi);
            T pmm0 = 1;
            T pmm1 = 1;
            for (index_t l_i = m; l_i < l; l_i++) {
                T P_lm = 1;
                
                // evaluate P_{l_i}_m(x) from previous values (if extant)
                if (l_i == m) {
                    if (m > 0) {
                        // compute double factorial of m, with (-1)^m factor
                        for (index_t i = 1; i <= m; i++) {
                            pmm0 *= -(i * 2 - 1);
                        }
                        pmm0 *= std::pow(1 - x * x, m / 2.0);
                    }
                    P_lm *= pmm0;
                } else if (l_i == m + 1) {
                    pmm1  = x * (2 * m + 1) * pmm0;
                    P_lm *= pmm1;
                } else {
                    T pmi = (x * (2 * l_i - 1) * pmm1 - (l_i + m - 1) * pmm0) / 
                            (l_i - m);
                    pmm0  = pmm1;
                    pmm1  = pmi;
                    P_lm *= pmi;
                }
                // P_lm now contains the legendre polynomial value at x
                
                // evaluate ((l-m)!/(l+m)!)^-1
                T z = 1;
                if (l_i != 0) for (index_t i  = l_i - m + 1; i <= l_i + m; i++) { z *= i; }
                
                // fold Legendre values into SH basis functions
                T k = std::sqrt((1/z) * (2 * l_i + 1) / (4 * pi));
                if (m != 0) {
                    sum += P_lm * sr2 * k * c_ma * coeff(l_i,  m);
                    sum += P_lm * sr2 * k * s_ma * coeff(l_i, -m);
                } else {
                    sum += P_lm * k * coeff(l_i, 0);
                }
            }
        }
        return sum;
    }
    
    
    void _project(T cos_alt, T azi, T val) {
        static const T sr2 = std::sqrt((T)2);
        static const T pi  = boost::math::constants::pi<T>();
        const index_t l = bands();
        const T x = cos_alt;
        for (index_t m = 0; m < l; m++) {
            T c_ma = std::cos(m * azi);
            T s_ma = std::sin(m * azi);
            T pmm0 = 1;
            T pmm1 = 1;
            for (index_t l_i = m; l_i < l; l_i++) {
                T P_lm = 1;
                // evaluate P_{l_i}_m(x) from previous values (if extant)
                if (l_i == m) {
                    if (m > 0) {
                        for (index_t i = 1; i <= m; i++) {
                            pmm0 *= -(i * 2 - 1);
                        }
                        pmm0 *= std::pow(1 - x * x, m / 2.0);
                    }
                    P_lm = pmm0;
                } else if (l_i == m + 1) {
                    pmm1  = x * (2 * m + 1) * pmm0;
                    P_lm = pmm1;
                } else {
                    T pmi = (x * (2 * l_i - 1) * pmm1 - (l_i + m - 1) * pmm0) / 
                            (l_i - m);
                    pmm0  = pmm1;
                    pmm1  = pmi;
                    P_lm  = pmi;
                }
                // P_lm now contains the legendre polynomial value at x
                
                // evaluate ((l-m)!/(l+m)!)^-1
                // this could be done with a LUT for less numerical precision
                // but I prefer "correct" over "moderately faster".
                T z = 1;
                if (l_i != 0) for (index_t i = l_i - m + 1; i <= l_i + m; i++) { z *= i; }
                
                // fold Legendre values into SH basis functions
                T k = std::sqrt((1/z) * (2 * l_i + 1) / (4 * pi));
                if (m != 0) {
                    coeff(l_i,  m) += P_lm * sr2 * k * c_ma * val;
                    coeff(l_i, -m) += P_lm * sr2 * k * s_ma * val;
                } else {
                    coeff(l_i,  0) += P_lm * k * val;
                }
            }
        }
    }
    
    public:
    
    /**
     * Evaluate this SphericalHarmonics function in given direction.
     * The polar direction is aligned with z+, and the x-axis has zero azimuthal 
     * angle.
     * 
     * @param d Direction along which to sample the funciton.
     * @return Reconstructed function value.
     */
    inline T eval(Vec<T,3> d) const {
        d = d.unit();
        return _eval(d.z, std::atan2(d.y, d.x));
    }
    
    /**
     * Evaluate this SphericalHarmonics function in the direction described by
     * `alt` and `azi`.
     * @param alt Angle from the polar axis.
     * @param azi Azimuthal (longitude) angle.
     * @return Reconstructed function value.
     */
    inline T eval(T alt, T azi) const {
        return _eval(std::cos(alt), azi);
    }
    
    
    /**
     * Project the Dirac delta function with direction `d` and integral `val` onto the 
     * spherical harmonic basis. The projection will be added to the current 
     * representation. In other words, "add" the supplied sample to the 
     * representation.
     * 
     * To represent a known function on the sphere, draw many samples from it 
     * and project each sample, then divide by 4&pi; times the number of samples.
     * 
     * @param alt Angle from the polar axis.
     * @param azi Azimuthal (longitude) angle.
     * @param val Value of sample.
     */
    inline void project(Vec<T,3> d, T val) {
        d = d.unit();
        _project(d.z, std::atan2(d.y, d.x), val);
    }
    
    /**
     * Project the Dirac delta function with direction given by `alt` and `azi`, 
     * and integral `val` onto the  spherical harmonic basis. The projection will 
     * be added to the current representation. In other words, "add" the supplied 
     * sample to the representation.
     * 
     * To represent a known function on the sphere, draw many samples from it 
     * and project each sample, then divide by the number of samples.
     * 
     * @param d Direction of sample to project. 
     * @param val Value of sample.
     */
    inline void project(T alt, T azi, T val) {
        _project(std::cos(alt), azi, val);
    }
    
    /**
     * Negate the values of this SphericalHarmonics function.
     */
    SphericalHarmonics<T,Bands>& operator-() {
        for (index_t i = 0; i < size(); i++) {
            coeffs[i] = -coeffs[i];
        }
        return *this;
    }
    
    /**
     * Multiply this SphericalHarmonics function by a constant value.
     * @param k Gain factor.
     */
    SphericalHarmonics<T,Bands>& operator*=(T k) {
        for (index_t i = 0; i < size(); i++) {
            coeffs[i] *= k;
        }
        return *this;
    }
    
    /**
     * Divide this SphericalHarmonics function by a constant value.
     * @param k Divisor.
     */
    SphericalHarmonics<T,Bands>& operator/=(T k) {
        for (index_t i = 0; i < size(); i++) {
            coeffs[i] /= k;
        }
        return *this;
    }
    
    /**
     * Add a constant value to this SphericalHarmonics function.
     * @param k Bias.
     */
    SphericalHarmonics<T,Bands>& operator+=(T k) {
        coeffs[0] += k;
        return *this;
    }
    
    /**
     * Add another SphericalHarmonics function to this one.
     */
    SphericalHarmonics<T,Bands>& operator+=(const SphericalHarmonics<T,Bands> &sh) {
        index_t n = std::min(sh.size(), size());
        for (index_t i = 0; i < n; i++) {
            coeffs[i] += sh.coeffs[i];
        }
        return *this;
    }
    
    /**
     * Construct a SphericalHarmonics function which is the sum of this and another.
     */
    SphericalHarmonics<T,Bands> operator+(const SphericalHarmonics<T,Bands> &sh) {
        index_t n = std::max(sh.size(), size());
        SphericalHarmonics<T,Bands> out(n);
        for (index_t i; i < n; i++) {
            out.coeffs[i] = coeffs[i] + sh.coeffs[i];
        }
    }
};

/// @} // group function
    
} // end namespace geom

#endif	/* SPHERICALHARMONICS_H */

