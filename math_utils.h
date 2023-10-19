#ifndef MATH_MATH_UTILS_H_
#define MATH_MATH_UTILS_H_

#include <limits.h>

/**
 * @brief Represents a small positive floating-point epsilon value, 1e-6 (1 ×
 * 10^(-6)).
 */
#define MU_EPS6 1e-6
/**
 * @brief Represents a small positive floating-point epsilon value, 1e-10 (1 ×
 * 10^(-10)).
 */
#define MU_EPS10 1e-10
/**
 * @brief Represents a small positive floating-point epsilon value, 1e-20 (1 ×
 * 10^(-20)).
 */
#define MU_EPS20 1e-20
/**
 * @brief Represents the exponential constant 1e10 (1 × 10^10).
 */
#define MU_E10 1e10
/**
 * @brief Represents the value of π (pi) accurate up to the specified number of
 * decimal places.
 */
#define MU_PI 3.14159265358979323846264338327950288
/**
 * @brief Represents the mathematical constant 'e', the base of natural
 * logarithms.
 */
#define MU_E 2.71828182845904523536028747135266250
/**
 * @brief Represents the natural logarithm of 2.
 */
#define MU_LN2 0.693147180559945309417232121458
/**
 * @brief Represents the logarithm of base 10.
 */
#define MU_LN10 2.30258509299404568401799145468
/**
 * @brief Represents the golden ratio, often denoted by the Greek letter phi
 * (φ).
 */
#define MU_PHI 1.6180339887498948482
/**
 * @brief Represents the reciprocal of the golden ratio (1/φ).
 */
#define MU_1_PHI 0.6180339887498948482
/**
 * @brief Represents the square root of 2.
 */
#define MU_SQRT2 1.4142135623730950488
/**
 * @brief Represents the square root of 3.
 */
#define MU_SQRT3 1.7320508075688772935
/**
 * @brief Represents the square root of 5.
 */
#define MU_SQRT5 2.23606797749978969
/**
 * @brief Represents the Catalan constant.
 */
#define MU_CATALAN 0.915965594177219
/**
 * @brief Represents Cahen's constant.
 */
#define MU_CAHEN 0.64341054629
/**
 * @brief Represents positive infinity.
 */
#define MU_INF 1.0 / 0.0
/**
 * @brief Represents a Not-a-Number (NaN) value.
 */
#define MU_NAN 0.0 / 0.0

/**
 * @brief Computes the absolute value of an integer.
 *
 * This function calculates the absolute value of the given integer.
 *
 * @param x Integer value for which absolute value is calculated.
 * @return Absolute value of the input integer.
 */
long int mu_abs(int x);

/**
 * @brief Computes the absolute value of a double-precision floating-point
 * number.
 *
 * This function calculates the absolute value of the given double-precision
 * floating-point number.
 *
 * @param x Double-precision floating-point number for which absolute value is
 * calculated.
 * @return Absolute value of the input double-precision floating-point number.
 */
long double mu_fabs(double x);

/**
 * @brief Truncates a double-precision floating-point number to an integer.
 *
 * This function truncates the given double-precision floating-point number `x`
 * to the nearest integer toward zero. It behaves like removing the fractional
 * part of the number, keeping only the integer portion.
 *
 * @param x Double-precision floating-point number to be truncated.
 * @return Nearest integer toward zero after truncating the input number.
 */
long double mu_trunc(double x);

/**
 * @brief Rounds up a double-precision floating-point number to the smallest
 * integral value not less than the input.
 *
 * This function rounds up the given double-precision floating-point number `x`
 * to the smallest integer that is not less than `x`. If `x` is already an
 * integer, it remains unchanged.
 *
 * @param x Double-precision floating-point number to be rounded up.
 * @return Smallest integer not less than the input `x`.
 */
long double mu_ceil(double x);

/**
 * @brief Rounds down a double-precision floating-point number to the largest
 * integral value not greater than the input.
 *
 * This function rounds down the given double-precision floating-point number
 * `x` to the largest integer that is not greater than `x`. If `x` is already an
 * integer, it remains unchanged.
 *
 * @param x Double-precision floating-point number to be rounded down.
 * @return Largest integer not greater than the input `x`.
 */
long double mu_floor(double x);

/**
 * @brief Computes the remainder of dividing two double-precision floating-point
 * numbers.
 *
 * This function calculates the remainder of dividing the absolute values of `x`
 * by `y`, preserving the sign of the dividend.
 *
 * @param x Double-precision floating-point dividend.
 * @param y Double-precision floating-point divisor.
 * @return Remainder of dividing `x` by `y`, with the same sign as `x`.
 */
long double mu_fmod(double x, double y);

/**
 * @brief Computes the sine of a double-precision floating-point number.
 *
 * This function calculates the sine of the given angle `x` in radians using the
 * Taylor series expansion.
 *
 * @param x Double-precision floating-point angle in radians.
 * @return Sine of the input angle `x`.
 */
long double mu_sin(double x);

/**
 * @brief Computes the cosine of a double-precision floating-point number.
 *
 * This function calculates the cosine of the given angle `x` in radians using
 * the Taylor series expansion.
 *
 * @param x Double-precision floating-point angle in radians.
 * @return Cosine of the input angle `x`.
 */
long double mu_cos(double x);

/**
 * @brief Computes the tangent of a double-precision floating-point number.
 *
 * This function calculates the tangent of the given angle `x` in radians by
 * the sine and the cosine Taylor series expansion.
 *
 * @param x Double-precision floating-point angle in radians.
 * @return Tangent of the input angle `x`.
 */
long double mu_tan(double x);

/**
 * @brief Computes the arcsine (inverse sine) of a double-precision
 * floating-point number.
 *
 * This function calculates the arcsine of the input `x`, returning the result
 * in radians using the Taylor series expansion.
 *
 * @param x Double-precision floating-point number for which the arcsine is
 * calculated.
 * @return Arcsine of the input `x` in radians.
 */
long double mu_asin(double x);

/**
 * @brief Computes the arccosine (inverse cosine) of a double-precision
 * floating-point number.
 *
 * This function calculates the arccosine of the input `x`, returning the result
 * in radians using the Taylor series expansion.
 *
 * @param x Double-precision floating-point number for which the arccosine is
 * calculated.
 * @return Arccosine of the input `x` in radians.
 */
long double mu_acos(double x);

/**
 * @brief Computes the arctangent (inverse tangent) of a double-precision
 * floating-point number.
 *
 * This function calculates the arctangent of the input `x`, returning the
 * result in radians using the Taylor series expansion.
 *
 * @param x Double-precision floating-point number for which the arctangent is
 * calculated.
 * @return Arctangent of the input `x` in radians.
 */
long double mu_atan(double x);

/**
 * @brief Computes the square root of a non-negative double-precision
 * floating-point number.
 *
 * This function calculates the square root of the input `x` using the Taylor
 * series expansion.
 *
 * @param x Non-negative double-precision floating-point number for which the
 * square root is calculated.
 * @return Square root of the input `x`.
 */
long double mu_sqrt(double x);

/**
 * @brief Computes the power of a double-precision floating-point number.
 *
 * This function calculates `base` raised to the power of `exp` using the Taylor
 * series expansion.
 *
 * @param base Base value.
 * @param exp Exponent value.
 * @return Result of `base` raised to the power of `exp`.
 */
long double mu_pow(double base, double exp);

/**
 * @brief Computes the exponential of a double-precision floating-point
 * number.
 *
 * This function calculates the exponential (e^x) of the given double-precision
 * floating-point number.
 *
 * @param x Double-precision floating-point number for which the exponential is
 * calculated.
 * @return Exponential of the input double-precision floating-point number.
 */
long double mu_exp(double x);

/**
 * @brief Computes the natural logarithm of a positive double-precision
 * floating-point number.
 *
 * This function calculates the natural (base e) logarithm of the given positive
 * double-precision floating-point number.
 *
 * @param x Positive double-precision floating-point number for which the
 * natural logarithm is calculated.
 * @return Natural logarithm of the input double-precision floating-point
 */
long double mu_log(double x);

#endif  // MATH_MATH_UTILS_H_
