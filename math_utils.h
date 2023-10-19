#ifndef MATH_MATH_UTILS_H_
#define MATH_MATH_UTILS_H_

#include <limits.h>

#define MU_EPS6 1e-6
#define MU_EPS10 1e-10
#define MU_EPS20 1e-20
#define MU_E10 1e10
#define MU_PI 3.14159265358979323846264338327950288
#define MU_E 2.71828182845904523536028747135266250
#define MU_LN2 0.693147180559945309417232121458
#define MU_LN10 2.30258509299404568401799145468
#define MU_GOLDEN 1.6180339887498948482
#define MU_REVGOLDEN 0.6180339887498948482
#define MU_SQRT2 1.4142135623730950488
#define MU_SQRT3 1.7320508075688772935
#define MU_SQRT5 2.23606797749978969
#define MU_CATALAN 0.915965594177219
#define MU_CAHEN 0.64341054629
#define MU_INF 1.0 / 0.0
#define MU_NAN 0.0 / 0.0

/**
 * @brief Calculates the absolute value of an integer.
 *
 * This function computes the absolute value of the given integer.
 *
 * @param x Integer value for which absolute value is calculated.
 * @return Absolute value of the input integer.
 */
long int mu_abs(int x);

/**
 * @brief Calculates the absolute value of a double-precision floating-point
 * number.
 *
 * This function computes the absolute value of the given double-precision
 * floating-point number.
 *
 * @param x Double-precision floating-point number for which absolute value is
 * calculated.
 * @return Absolute value of the input double-precision floating-point number.
 */
long double mu_fabs(double x);
long double mu_trunc(double x);
long double mu_ceil(double x);
long double mu_floor(double x);
long double mu_fmod(double x, double y);
long double mu_sin(double x);
long double mu_cos(double x);
long double mu_tan(double x);
long double mu_asin(double x);

/**
 * @brief Calculates the arccosine (inverse cosine) of a given double value.
 *
 * This function computes the arccosine of the input, returning the result in
 * radians.
 *
 * @param x Double value for which the arccosine is calculated.
 * @return Arccosine of the input in radians. If is outside the range [-1, 1],
 * the function returns NaN (Not a Number).
 */
long double mu_acos(double x);
long double mu_atan(double x);
long double mu_sqrt(double x);
long double mu_pow(double base, double exp);
long double mu_exp(double x);
long double mu_log(double x);

#endif  // MATH_MATH_UTILS_H_
