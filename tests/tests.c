#include <check.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "math_utils.h"

void run_range_tests(long double (*mu_func)(double), double (*std_func)(double),
                     double start, double end, double step, double tolerance) {
  for (double x = start; x < end; x += step) {
    ck_assert_ldouble_eq_tol(mu_func(x), std_func(x), tolerance);
  }
}

void run_const_tests(long double (*mu_func)(double), double (*std_func)(double),
                     double tolerance) {
  double constants[] = {MU_PI,    MU_LN10,  MU_E,       MU_SQRT2,
                        MU_SQRT3, MU_SQRT5, MU_CATALAN, MU_CAHEN,
                        MU_LN2,   MU_PHI,   MU_1_PHI};

  for (size_t i = 0; i < sizeof(constants) / sizeof(constants[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_func(constants[i]), std_func(constants[i]),
                             tolerance);
  }
}

void run_random_tests(long double (*mu_func)(double),
                      double (*std_func)(double), double low, double high,
                      double tolerance) {
  srand((unsigned int)time(NULL));
  for (int i = 0; i < 1000; ++i) {
    double x = low + (double)rand() / RAND_MAX * (high - low);
    ck_assert_ldouble_eq_tol(mu_func(x), std_func(x), tolerance);
  }
}

void run_const_tests_2args(long double (*mu_func)(double, double),
                           double (*std_func)(double, double),
                           double tolerance) {
  double constants[] = {MU_PI,    MU_LN10,  MU_E,       MU_SQRT2,
                        MU_SQRT3, MU_SQRT5, MU_CATALAN, MU_CAHEN,
                        MU_LN2,   MU_PHI,   MU_1_PHI};

  for (size_t i = 0; i < sizeof(constants) / sizeof(constants[0]); ++i) {
    for (size_t j = 0; j < sizeof(constants) / sizeof(constants[0]); ++j) {
      ck_assert_ldouble_eq_tol(mu_func(constants[i], constants[j]),
                               std_func(constants[i], constants[j]), tolerance);
    }
  }
}

void run_random_tests_2args(long double (*mu_func)(double, double),
                            double (*std_func)(double, double), double low1,
                            double high1, double low2, double high2,
                            double tolerance) {
  srand((unsigned int)time(NULL));
  for (int i = 0; i < 1000; ++i) {
    double x = low1 + (double)rand() / RAND_MAX * (high1 - low1);
    double y = low2 + (double)rand() / RAND_MAX * (high2 - low2);
    ck_assert_ldouble_eq_tol(mu_func(x, y), std_func(x, y), tolerance);
  }
}

START_TEST(test_mu_abs) {
  for (int x = -10000; x < 10000; ++x) {
    ck_assert_int_eq(mu_abs(x), abs(x));
  }

  for (unsigned int i = 0; i < 1000; ++i) {
    int x = rand_r(&i);
    ck_assert_int_eq(mu_abs(x), abs(x));
  }
}
END_TEST

START_TEST(test_mu_fabs) {
  run_range_tests(mu_fabs, fabs, -1000.0, 1000.0, 0.1, MU_EPS6);
  run_const_tests(mu_fabs, fabs, MU_EPS6);
  run_random_tests(mu_fabs, fabs, -MU_E10, MU_E10, MU_EPS6);

  ck_assert_ldouble_nan(mu_fabs(MU_NAN));
  ck_assert_ldouble_eq(mu_fabs(MU_INF), fabs(MU_INF));
  ck_assert_ldouble_eq(mu_fabs(-MU_INF), fabs(-MU_INF));
}
END_TEST

START_TEST(test_mu_trunc) {
  run_range_tests(mu_trunc, trunc, -1000.0, 1000.0, 0.1, MU_EPS6);
  run_const_tests(mu_trunc, trunc, MU_EPS6);
  run_random_tests(mu_trunc, trunc, -MU_E10, MU_E10, MU_EPS6);

  ck_assert_ldouble_nan(mu_trunc(MU_NAN));
  ck_assert_ldouble_eq(mu_trunc(MU_INF), trunc(MU_INF));
  ck_assert_ldouble_eq(mu_trunc(-MU_INF), trunc(-MU_INF));
}
END_TEST

START_TEST(test_mu_ceil) {
  run_range_tests(mu_ceil, ceil, -1000.0, 1000.0, 0.1, MU_EPS6);
  run_const_tests(mu_ceil, ceil, MU_EPS6);
  run_random_tests(mu_ceil, ceil, -MU_E10, MU_E10, MU_EPS6);

  ck_assert_ldouble_nan(mu_ceil(MU_NAN));
  ck_assert_ldouble_eq(mu_ceil(MU_INF), ceil(MU_INF));
  ck_assert_ldouble_eq(mu_ceil(-MU_INF), ceil(-MU_INF));
}
END_TEST

START_TEST(test_mu_floor) {
  run_range_tests(mu_floor, floor, -1000.0, 1000.0, 0.1, MU_EPS6);
  run_const_tests(mu_floor, floor, MU_EPS6);
  run_random_tests(mu_floor, floor, -MU_E10, MU_E10, MU_EPS6);

  ck_assert_ldouble_nan(mu_floor(MU_NAN));
  ck_assert_ldouble_eq(mu_floor(MU_INF), floor(MU_INF));
  ck_assert_ldouble_eq(mu_floor(-MU_INF), floor(-MU_INF));
}
END_TEST

START_TEST(test_mu_fmod) {
  run_const_tests_2args(mu_fmod, fmod, MU_EPS6);
  run_random_tests_2args(mu_fmod, fmod, -MU_E10, MU_E10, -MU_E10, MU_E10,
                         MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_fmod(MU_PI, MU_E), fmod(MU_PI, MU_E), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_fmod(MU_E, MU_PI), fmod(MU_E, MU_PI), MU_EPS6);

  ck_assert_ldouble_nan(mu_fmod(MU_INF, MU_INF));
  ck_assert_ldouble_nan(mu_fmod(MU_INF, 0.0));
  ck_assert_ldouble_nan(mu_fmod(MU_INF, 1.0));
  ck_assert_ldouble_nan(mu_fmod(MU_NAN, MU_NAN));
  ck_assert_ldouble_nan(mu_fmod(MU_NAN, 1.0));
  ck_assert_ldouble_nan(mu_fmod(1.0, MU_NAN));
  ck_assert_ldouble_nan(mu_fmod(1.0, 0.0));
  ck_assert_ldouble_nan(mu_fmod(0.0, 0.0));
  ck_assert_ldouble_eq(mu_fmod(1.0, MU_INF), fmod(1.0, MU_INF));
  ck_assert_ldouble_eq(mu_fmod(MU_E10, -MU_INF), fmod(MU_E10, -MU_INF));
}
END_TEST

START_TEST(test_mu_sin) {
  run_range_tests(mu_sin, sin, -1000.0, 1000.0, 0.1, MU_EPS6);
  run_const_tests(mu_sin, sin, MU_EPS6);
  run_random_tests(mu_sin, sin, -MU_E10, MU_E10, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_sin(0.0), sin(0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(MU_PI / 2), sin(MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(-MU_PI / 2), sin(-MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(MU_PI), sin(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(-MU_PI), sin(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(2 * MU_PI), sin(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(-2 * MU_PI), sin(-2 * MU_PI), MU_EPS6);

  ck_assert_ldouble_nan(mu_sin(MU_NAN));
  ck_assert_ldouble_nan(mu_sin(MU_INF));
  ck_assert_ldouble_nan(mu_sin(-MU_INF));
}
END_TEST

START_TEST(test_mu_cos) {
  run_range_tests(mu_cos, cos, -1000.0, 1000.0, 0.1, MU_EPS6);
  run_const_tests(mu_cos, cos, MU_EPS6);
  run_random_tests(mu_cos, cos, -MU_E10, MU_E10, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_cos(0), cos(0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(MU_PI / 2), cos(MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(-MU_PI / 2), cos(-MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(MU_PI), cos(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(-MU_PI), cos(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(2 * MU_PI), cos(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(-2 * MU_PI), cos(-2 * MU_PI), MU_EPS6);

  ck_assert_ldouble_nan(mu_cos(MU_NAN));
  ck_assert_ldouble_nan(mu_cos(MU_INF));
  ck_assert_ldouble_nan(mu_cos(-MU_INF));
}
END_TEST

START_TEST(test_mu_tan) {
  run_range_tests(mu_tan, tan, -1.0, 1.0, 0.002, MU_EPS6);
  run_const_tests(mu_tan, tan, MU_EPS6);
  run_random_tests(mu_cos, cos, -1.0, 1.0, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_tan(0), tan(0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(MU_PI), tan(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(-MU_PI), tan(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(2 * MU_PI), tan(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(-2 * MU_PI), tan(-2 * MU_PI), MU_EPS6);

  ck_assert_ldouble_nan(mu_tan(MU_NAN));
  ck_assert_ldouble_nan(mu_tan(MU_INF));
  ck_assert_ldouble_nan(mu_tan(-MU_INF));
}
END_TEST

START_TEST(test_mu_asin) {
  run_range_tests(mu_asin, asin, -1.0, 1.0, 0.002, MU_EPS6);
  run_random_tests(mu_asin, asin, -0.999, 0.999, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_asin(0.0), asin(0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_asin(MU_LN2), asin(MU_LN2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_asin(MU_1_PHI), asin(MU_1_PHI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_asin(MU_CATALAN), asin(MU_CATALAN), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_asin(MU_CAHEN), asin(MU_CAHEN), MU_EPS6);

  ck_assert_ldouble_nan(mu_asin(MU_PI));
  ck_assert_ldouble_nan(mu_asin(MU_E));
  ck_assert_ldouble_nan(mu_asin(MU_LN10));
  ck_assert_ldouble_nan(mu_asin(MU_PHI));
  ck_assert_ldouble_nan(mu_asin(MU_SQRT2));
  ck_assert_ldouble_nan(mu_asin(MU_SQRT3));
  ck_assert_ldouble_nan(mu_asin(MU_SQRT5));

  ck_assert_ldouble_nan(mu_asin(MU_NAN));
  ck_assert_ldouble_nan(mu_asin(MU_INF));
  ck_assert_ldouble_nan(mu_asin(-MU_INF));
}
END_TEST

START_TEST(test_mu_acos) {
  run_range_tests(mu_acos, acos, -1.0, 1.0, 0.002, MU_EPS6);
  run_random_tests(mu_acos, acos, -0.999, 0.999, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_acos(0.0), acos(0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_acos(0.0), acos(0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_acos(MU_LN2), acos(MU_LN2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_acos(MU_1_PHI), acos(MU_1_PHI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_acos(MU_CATALAN), acos(MU_CATALAN), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_acos(MU_CAHEN), acos(MU_CAHEN), MU_EPS6);

  ck_assert_ldouble_nan(mu_acos(MU_PI));
  ck_assert_ldouble_nan(mu_acos(MU_E));
  ck_assert_ldouble_nan(mu_acos(MU_LN10));
  ck_assert_ldouble_nan(mu_acos(MU_PHI));
  ck_assert_ldouble_nan(mu_acos(MU_SQRT2));
  ck_assert_ldouble_nan(mu_acos(MU_SQRT3));
  ck_assert_ldouble_nan(mu_acos(MU_SQRT5));

  ck_assert_ldouble_nan(mu_acos(MU_NAN));
  ck_assert_ldouble_nan(mu_acos(MU_INF));
  ck_assert_ldouble_nan(mu_acos(-MU_INF));
}
END_TEST

START_TEST(test_mu_atan) {
  run_range_tests(mu_atan, atan, -10.0, 10.0, 0.1, MU_EPS6);
  run_const_tests(mu_atan, atan, MU_EPS6);
  run_random_tests(mu_atan, atan, -10.0, 10.0, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_atan(0.0), atan(0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(MU_PI / 2), atan(MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-MU_PI / 2), atan(-MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(MU_PI), atan(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-MU_PI), atan(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(2 * MU_PI), atan(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-2 * MU_PI), atan(-2 * MU_PI), MU_EPS6);

  ck_assert_ldouble_nan(mu_atan(MU_NAN));
  ck_assert_ldouble_eq_tol(mu_atan(MU_INF), atan(MU_INF), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-MU_INF), atan(-MU_INF), MU_EPS6);
}
END_TEST

START_TEST(test_mu_sqrt) {
  run_range_tests(mu_sqrt, sqrt, 0.0, 10000.0, 10, MU_EPS6);
  run_range_tests(mu_sqrt, sqrt, 0.0, 1.0, 0.001, MU_EPS6);
  run_const_tests(mu_sqrt, sqrt, MU_EPS6);
  run_random_tests(mu_sqrt, sqrt, 0.0, MU_E10, MU_EPS6);

  ck_assert_ldouble_eq_tol(mu_sqrt(-0.0), sqrt(-0.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_sqrt(-MU_EPS6));
  ck_assert_ldouble_nan(mu_sqrt(MU_NAN));
  ck_assert_ldouble_nan(mu_sqrt(MU_INF));
  ck_assert_ldouble_nan(mu_sqrt(-MU_INF));
}
END_TEST

START_TEST(test_mu_pow) {
  run_const_tests_2args(mu_pow, pow, MU_EPS6);
  run_random_tests_2args(mu_pow, pow, 1, 40, -5, 5, MU_EPS6);

  ck_assert_ldouble_eq(mu_pow(-MU_INF, -MU_INF), pow(-MU_INF, -MU_INF));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, -10.0), pow(-MU_INF, -10.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, -1.0), pow(-MU_INF, -1.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, -0.1), pow(-MU_INF, -0.1));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 0.0), pow(-MU_INF, 0.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 0.1), pow(-MU_INF, 0.1));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 1.0), pow(-MU_INF, 1.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 9.0), pow(-MU_INF, 9.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 10.0), pow(-MU_INF, 10.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, MU_INF), pow(-MU_INF, MU_INF));
  ck_assert_ldouble_nan(mu_pow(-MU_INF, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(-10.0, -MU_INF), pow(-10.0, -MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(-10.0, -10.0), pow(-10.0, -10.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-10.0, -1.0), pow(-10.0, -1.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(-10.0, -0.1));
  ck_assert_ldouble_eq_tol(mu_pow(-10.0, 0.0), pow(-10.0, 0.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(-10.0, 0.1));
  ck_assert_ldouble_eq_tol(mu_pow(-10.0, 1.0), pow(-10.0, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-10.0, 9.0), pow(-10.0, 9.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-10.0, 10.0), pow(-10.0, 10.0), MU_EPS6);
  ck_assert_ldouble_eq(mu_pow(-10.0, MU_INF), pow(-10.0, MU_INF));
  ck_assert_ldouble_nan(mu_pow(-10.0, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(-1.0, -MU_INF), pow(-1.0, -MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(-1.0, -10.0), pow(-1.0, -10.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-1.0, -1.0), pow(-1.0, -1.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(-1.0, -0.1));
  ck_assert_ldouble_eq_tol(mu_pow(-1.0, 0.0), pow(-1.0, 0.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(-1.0, 0.1));
  ck_assert_ldouble_eq_tol(mu_pow(-1.0, 1.0), pow(-1.0, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-1.0, 9.0), pow(-1.0, 9.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-1.0, 10.0), pow(-1.0, 10.0), MU_EPS6);
  ck_assert_ldouble_eq(mu_pow(-1.0, MU_INF), pow(-1.0, MU_INF));
  ck_assert_ldouble_nan(mu_pow(-1.0, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(-0.1, -MU_INF), pow(-0.1, -MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(-0.1, -10.0), pow(-0.1, -10.0), MU_E10);
  ck_assert_ldouble_eq_tol(mu_pow(-0.1, -1.0), pow(-0.1, -1.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(-0.1, -0.1));
  ck_assert_ldouble_eq_tol(mu_pow(-0.1, 0.0), pow(-0.1, 0.0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(-0.1, 0.1));
  ck_assert_ldouble_eq_tol(mu_pow(-0.1, 1.0), pow(-0.1, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(-0.1, 9.0), pow(-0.1, 9.0), MU_EPS10);
  ck_assert_ldouble_eq_tol(mu_pow(-0.1, 10.0), pow(-0.1, 10.0), MU_EPS10);
  ck_assert_ldouble_eq(mu_pow(-0.1, MU_INF), pow(-0.1, MU_INF));
  ck_assert_ldouble_nan(mu_pow(-0.1, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(0.0, -MU_INF), pow(0.0, -MU_INF));
  ck_assert_ldouble_eq(mu_pow(0.0, -10.0), pow(0.0, -10.0));
  ck_assert_ldouble_eq(mu_pow(0.0, -1.0), pow(0.0, -1.0));
  ck_assert_ldouble_eq(mu_pow(0.0, -0.1), pow(0.0, -0.1));
  ck_assert_ldouble_eq_tol(mu_pow(0.0, 0.0), pow(0.0, 0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.0, 0.1), pow(0.0, 0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.0, 1.0), pow(0.0, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.0, 10.0), pow(0.0, 10.0), MU_EPS6);
  ck_assert_ldouble_eq(mu_pow(0.0, MU_INF), pow(0.0, MU_INF));
  ck_assert_ldouble_nan(mu_pow(0.0, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(0.1, -MU_INF), pow(0.1, -MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(0.1, -10.0), pow(0.1, -10.0), MU_E10);
  ck_assert_ldouble_eq_tol(mu_pow(0.1, -1.0), pow(0.1, -1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.1, -0.1), pow(0.1, -0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.1, 0.0), pow(0.1, 0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.1, 0.1), pow(0.1, 0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.1, 1.0), pow(0.1, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(0.1, 10.0), pow(0.1, 10.0), MU_EPS10);
  ck_assert_ldouble_eq(mu_pow(0.1, MU_INF), pow(0.1, MU_INF));
  ck_assert_ldouble_nan(mu_pow(0.1, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(1.0, -MU_INF), pow(1.0, -MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(1.0, -10.0), pow(1.0, -10.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(1.0, -1.0), pow(1.0, -1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(1.0, -0.1), pow(1.0, -0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(1.0, 0.0), pow(1.0, 0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(1.0, 0.1), pow(1.0, 0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(1.0, 1.0), pow(1.0, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(1.0, 10.0), pow(1.0, 10.0), MU_EPS6);
  ck_assert_ldouble_eq(mu_pow(1.0, MU_INF), pow(1.0, MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(1.0, MU_NAN), pow(1.0, MU_NAN), MU_EPS6);

  ck_assert_ldouble_eq(mu_pow(10.0, -MU_INF), pow(10.0, -MU_INF));
  ck_assert_ldouble_eq_tol(mu_pow(10.0, -10.0), pow(10.0, -10.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(10.0, -1.0), pow(10.0, -1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(10.0, -0.1), pow(10.0, -0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(10.0, 0.0), pow(10.0, 0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(10.0, 0.1), pow(10.0, 0.1), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(10.0, 1.0), pow(10.0, 1.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_pow(10.0, 10.0), pow(10.0, 10.0), 1);
  ck_assert_ldouble_eq(mu_pow(10.0, MU_INF), pow(10.0, MU_INF));
  ck_assert_ldouble_nan(mu_pow(10.0, MU_NAN));

  ck_assert_ldouble_eq(mu_pow(MU_INF, -MU_INF), pow(MU_INF, -MU_INF));
  ck_assert_ldouble_eq(mu_pow(MU_INF, -10.0), pow(MU_INF, -10.0));
  ck_assert_ldouble_eq(mu_pow(MU_INF, -1.0), pow(MU_INF, -1.0));
  ck_assert_ldouble_eq(mu_pow(MU_INF, -0.1), pow(MU_INF, -0.1));
  ck_assert_ldouble_eq(mu_pow(MU_INF, 0.0), pow(MU_INF, 0.0));
  ck_assert_ldouble_eq(mu_pow(MU_INF, 0.1), pow(MU_INF, 0.1));
  ck_assert_ldouble_eq(mu_pow(MU_INF, 1.0), pow(MU_INF, 1.0));
  ck_assert_ldouble_eq(mu_pow(MU_INF, 10.0), pow(MU_INF, 10.0));
  ck_assert_ldouble_eq(mu_pow(MU_INF, MU_INF), pow(MU_INF, MU_INF));
  ck_assert_ldouble_nan(mu_pow(MU_INF, MU_NAN));

  ck_assert_ldouble_nan(mu_pow(MU_NAN, -MU_INF));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, -10.0));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, -1.0));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, -0.1));
  ck_assert_ldouble_eq_tol(mu_pow(MU_NAN, 0), pow(MU_NAN, 0), MU_EPS6);
  ck_assert_ldouble_nan(mu_pow(MU_NAN, 0.1));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, 1.0));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, 10.0));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, MU_INF));
  ck_assert_ldouble_nan(mu_pow(MU_NAN, MU_NAN));
}
END_TEST

START_TEST(test_mu_exp) {
  run_range_tests(mu_exp, exp, -100.0, 20.0, 0.1, MU_EPS6);
  run_const_tests(mu_exp, exp, MU_EPS6);
  run_random_tests(mu_exp, exp, -50.0, 20, MU_EPS6);

  ck_assert_ldouble_nan(mu_exp(MU_NAN));
  ck_assert_ldouble_eq(mu_exp(MU_INF), exp(MU_INF));
  ck_assert_ldouble_eq_tol(mu_exp(-MU_INF), exp(-MU_INF), MU_EPS6);
}
END_TEST

START_TEST(test_mu_log) {
  run_range_tests(mu_log, log, 0.01, 2.0, 0.01, MU_EPS6);
  run_const_tests(mu_log, log, MU_EPS6);
  run_random_tests(mu_log, log, MU_EPS20, MU_E10, MU_EPS6);

  ck_assert_ldouble_nan(mu_log(-MU_EPS6));
  ck_assert_ldouble_nan(mu_log(MU_NAN));
  ck_assert_ldouble_nan(mu_log(-MU_INF));
  ck_assert_ldouble_eq(mu_log(0.0), log(0.0));
  ck_assert_ldouble_eq(mu_log(MU_INF), log(MU_INF));
}
END_TEST

Suite *math_utils_suite(void) {
  Suite *suite;
  TCase *core;

  suite = suite_create("math_utils");
  core = tcase_create("Core");

  tcase_add_test(core, test_mu_abs);
  tcase_add_test(core, test_mu_fabs);
  tcase_add_test(core, test_mu_trunc);
  tcase_add_test(core, test_mu_ceil);
  tcase_add_test(core, test_mu_floor);
  tcase_add_test(core, test_mu_fmod);
  tcase_add_test(core, test_mu_sin);
  tcase_add_test(core, test_mu_cos);
  tcase_add_test(core, test_mu_tan);
  tcase_add_test(core, test_mu_asin);
  tcase_add_test(core, test_mu_acos);
  tcase_add_test(core, test_mu_atan);
  tcase_add_test(core, test_mu_sqrt);
  tcase_add_test(core, test_mu_pow);
  tcase_add_test(core, test_mu_exp);
  tcase_add_test(core, test_mu_log);

  suite_add_tcase(suite, core);

  return (suite);
}

int main(void) {
  int failed = 0;
  Suite *suite;

  SRunner *runner;

  suite = math_utils_suite();
  runner = srunner_create(suite);

  srunner_run_all(runner, CK_NORMAL);
  failed = srunner_ntests_failed(runner);
  srunner_free(runner);

  return (failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE);
}
