#include <check.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "math_utils.h"

void init_random() { srand((unsigned int)time(NULL)); }

double rand_double(double low, double high) {
  return low + (double)rand() / RAND_MAX * (high - low);
}

double test_cases[15] = {MU_PI,    MU_EPS20, MU_LN10,    MU_E,     MU_SQRT2,
                         MU_SQRT3, MU_SQRT5, MU_CATALAN, MU_CAHEN, MU_LN2,
                         MU_PHI,   MU_1_PHI, 0.0,        -0.0};

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
  for (double x = -1000.0; x < 1000.0; x += 0.1) {
    ck_assert_ldouble_eq(mu_fabs(x), fabs(x));
  }

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq(mu_fabs(test_cases[i]), fabs(test_cases[i]));
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, MU_E10);
    ck_assert_ldouble_eq(mu_fabs(x), fabs(x));
  }

  ck_assert_ldouble_nan(mu_fabs(MU_NAN));
  ck_assert_ldouble_eq(mu_fabs(MU_INF), fabs(MU_INF));
  ck_assert_ldouble_eq(mu_fabs(-MU_INF), fabs(-MU_INF));
}
END_TEST

START_TEST(test_mu_trunc) {
  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq(mu_trunc(test_cases[i]), trunc(test_cases[i]));
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, MU_E10);
    ck_assert_ldouble_eq(mu_trunc(x), trunc(x));
  }

  ck_assert_ldouble_nan(mu_trunc(MU_NAN));
  ck_assert_ldouble_eq(mu_trunc(MU_INF), trunc(MU_INF));
  ck_assert_ldouble_eq(mu_trunc(-MU_INF), trunc(-MU_INF));
}
END_TEST

START_TEST(test_mu_ceil) {
  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq(mu_ceil(test_cases[i]), ceil(test_cases[i]));
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, MU_E10);
    ck_assert_ldouble_eq(mu_ceil(x), ceil(x));
  }

  ck_assert_ldouble_nan(mu_ceil(MU_NAN));
  ck_assert_ldouble_eq(mu_ceil(MU_INF), ceil(MU_INF));
  ck_assert_ldouble_eq(mu_ceil(-MU_INF), ceil(-MU_INF));
}
END_TEST

START_TEST(test_mu_floor) {
  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq(mu_floor(test_cases[i]), floor(test_cases[i]));
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, MU_E10);
    ck_assert_ldouble_eq(mu_floor(x), floor(x));
  }

  ck_assert_ldouble_nan(mu_floor(MU_NAN));
  ck_assert_ldouble_eq(mu_floor(MU_INF), floor(MU_INF));
  ck_assert_ldouble_eq(mu_floor(-MU_INF), floor(-MU_INF));
}
END_TEST

START_TEST(test_mu_fmod) {
  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, MU_E10);
    double y = rand_double(-MU_E10, MU_E10);
    ck_assert_ldouble_eq_tol(mu_fmod(x, y), fmod(x, y), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_fmod(MU_PI, MU_E), fmod(MU_PI, MU_E), MU_EPS6);

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
  for (double x = -1000.0; x < 1000.0; x += 0.1) {
    ck_assert_ldouble_eq_tol(mu_sin(x), sin(x), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_sin(0), sin(0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(MU_PI / 2), sin(MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(-MU_PI / 2), sin(-MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(MU_PI), sin(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(-MU_PI), sin(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(2 * MU_PI), sin(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_sin(-2 * MU_PI), sin(-2 * MU_PI), MU_EPS6);

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_sin(test_cases[i]), sin(test_cases[i]),
                             MU_EPS6);
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-2 * MU_PI, 2 * MU_PI);
    ck_assert_ldouble_eq_tol(mu_sin(x), sin(x), MU_EPS6);
  }

  ck_assert_ldouble_nan(mu_sin(MU_NAN));
  ck_assert_ldouble_nan(mu_sin(MU_INF));
  ck_assert_ldouble_nan(mu_sin(-MU_INF));
}
END_TEST

START_TEST(test_mu_cos) {
  for (double x = -1000.0; x < 1000.0; x += 0.1) {
    ck_assert_ldouble_eq_tol(mu_cos(x), cos(x), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_cos(0), cos(0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(MU_PI / 2), cos(MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(-MU_PI / 2), cos(-MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(MU_PI), cos(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(-MU_PI), cos(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(2 * MU_PI), cos(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_cos(-2 * MU_PI), cos(-2 * MU_PI), MU_EPS6);

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_cos(test_cases[i]), cos(test_cases[i]),
                             MU_EPS6);
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-2 * MU_PI, 2 * MU_PI);
    ck_assert_ldouble_eq_tol(mu_cos(x), cos(x), MU_EPS6);
  }

  ck_assert_ldouble_nan(mu_cos(MU_NAN));
  ck_assert_ldouble_nan(mu_cos(MU_INF));
  ck_assert_ldouble_nan(mu_cos(-MU_INF));
}
END_TEST

START_TEST(test_mu_tan) {
  for (double x = -1.0; x < 1.0; x += 0.002) {
    ck_assert_ldouble_eq_tol(mu_tan(x), tan(x), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_tan(0), tan(0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(MU_PI), tan(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(-MU_PI), tan(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(2 * MU_PI), tan(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_tan(-2 * MU_PI), tan(-2 * MU_PI), MU_EPS6);

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_tan(test_cases[i]), tan(test_cases[i]),
                             MU_EPS6);
  }

  ck_assert_ldouble_nan(mu_tan(MU_NAN));
  ck_assert_ldouble_nan(mu_tan(MU_INF));
  ck_assert_ldouble_nan(mu_tan(-MU_INF));
}
END_TEST

START_TEST(test_mu_asin) {
  for (double x = -1.0; x <= 1.0; x += 0.002) {
    ck_assert_ldouble_eq_tol(mu_asin(x), asin(x), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_asin(0.0), asin(0.0), MU_EPS6);

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    if (test_cases[i] >= -1 && test_cases[i] <= 1) {
      ck_assert_ldouble_eq_tol(mu_asin(test_cases[i]), asin(test_cases[i]),
                               MU_EPS6);
    } else {
      ck_assert_ldouble_nan(mu_asin(test_cases[i]));
    }
  }

  ck_assert_ldouble_nan(mu_asin(MU_NAN));
  ck_assert_ldouble_nan(mu_asin(MU_INF));
  ck_assert_ldouble_nan(mu_asin(-MU_INF));
}
END_TEST

START_TEST(test_mu_acos) {
  for (double x = -1.0; x <= 1.0; x += 0.002) {
    ck_assert_ldouble_eq_tol(mu_acos(x), acos(x), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_acos(0.0), acos(0.0), MU_EPS6);

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    if (test_cases[i] >= -1 && test_cases[i] <= 1) {
      ck_assert_ldouble_eq_tol(mu_acos(test_cases[i]), acos(test_cases[i]),
                               MU_EPS6);
    } else {
      ck_assert_ldouble_nan(mu_acos(test_cases[i]));
    }
  }

  ck_assert_ldouble_nan(mu_acos(MU_NAN));
  ck_assert_ldouble_nan(mu_acos(MU_INF));
  ck_assert_ldouble_nan(mu_acos(-MU_INF));
}
END_TEST

START_TEST(test_mu_atan) {
  for (double x = -10.0; x < 10.0; x += 0.1) {
    ck_assert_ldouble_eq_tol(mu_atan(x), atan(x), MU_EPS6);
  }

  ck_assert_ldouble_eq_tol(mu_atan(0.0), atan(0.0), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(MU_PI / 2), atan(MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-MU_PI / 2), atan(-MU_PI / 2), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(MU_PI), atan(MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-MU_PI), atan(-MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(2 * MU_PI), atan(2 * MU_PI), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-2 * MU_PI), atan(-2 * MU_PI), MU_EPS6);

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_atan(test_cases[i]), atan(test_cases[i]),
                             MU_EPS6);
  }

  ck_assert_ldouble_nan(mu_atan(MU_NAN));
  ck_assert_ldouble_eq_tol(mu_atan(MU_INF), atan(MU_INF), MU_EPS6);
  ck_assert_ldouble_eq_tol(mu_atan(-MU_INF), atan(-MU_INF), MU_EPS6);
}
END_TEST

START_TEST(test_mu_sqrt) {
  for (double x = -10000.0; x < 10000.0; x += 10) {
    if (x >= 0) {
      ck_assert_ldouble_eq_tol(mu_sqrt(x), sqrt(x), MU_EPS6);
    } else {
      ck_assert_ldouble_nan(mu_sqrt(x));
    }
  }

  for (double i = -1.0; i < 1.0; i += 0.001) {
    if (i >= 0) {
      ck_assert_ldouble_eq_tol(mu_sqrt(i), sqrt(i), MU_EPS6);
    } else {
      ck_assert_ldouble_nan(mu_sqrt(i));
    }
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, MU_E10);
    if (x >= 0)
      ck_assert_ldouble_eq_tol(mu_sqrt(x), sqrt(x), MU_EPS6);
    else
      ck_assert_ldouble_nan(mu_sqrt(x));
  }

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_sqrt(test_cases[i]), sqrt(test_cases[i]),
                             MU_EPS6);
  }

  ck_assert_ldouble_nan(mu_sqrt(MU_NAN));
  ck_assert_ldouble_nan(mu_sqrt(MU_INF));
  ck_assert_ldouble_nan(mu_sqrt(-MU_INF));
}
END_TEST

START_TEST(test_mu_pow) {
  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(1, 40);
    double y = rand_double(-5, 5);
    ck_assert_ldouble_eq_tol(mu_pow(x, y), pow(x, y), MU_EPS6);
  }

  for (int i = 0; i < 10; ++i) {
    double x = rand_double(-10, -1);
    ck_assert_ldouble_eq_tol(mu_pow(x, i), pow(x, i), MU_EPS6);
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-40, -1);
    double y = rand_double(-5, 5);
    ck_assert_ldouble_nan(mu_pow(x, y));
  }

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_pow(test_cases[i], i), pow(test_cases[i], i),
                             MU_EPS6);
    ck_assert_ldouble_eq_tol(mu_pow(i, test_cases[i]), pow(i, test_cases[i]),
                             MU_EPS6);
    if (test_cases[i] > MU_EPS6) {
      ck_assert_ldouble_eq_tol(mu_pow(test_cases[i], test_cases[i]),
                               pow(test_cases[i], test_cases[i]), MU_EPS6);
    }
  }

  ck_assert_ldouble_eq(mu_pow(-MU_INF, -MU_INF), pow(-MU_INF, -MU_INF));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, -10.0), pow(-MU_INF, -10.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, -1.0), pow(-MU_INF, -1.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, -0.1), pow(-MU_INF, -0.1));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 0.0), pow(-MU_INF, 0.0));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 0.1), pow(-MU_INF, 0.1));
  ck_assert_ldouble_eq(mu_pow(-MU_INF, 1.0), pow(-MU_INF, 1.0));
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
  for (double x = -100.0; x < 20.0; x += 0.1) {
    ck_assert_ldouble_eq_tol(mu_exp(x), exp(x), MU_EPS6);
  }

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    ck_assert_ldouble_eq_tol(mu_exp(test_cases[i]), exp(test_cases[i]),
                             MU_EPS6);
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-50, 20);
    ck_assert_ldouble_eq_tol(mu_exp(x), exp(x), MU_EPS6);
  }

  ck_assert_ldouble_nan(mu_exp(MU_NAN));
  ck_assert_ldouble_eq(mu_exp(MU_INF), exp(MU_INF));
  ck_assert_ldouble_eq_tol(mu_exp(-MU_INF), exp(-MU_INF), MU_EPS6);
}
END_TEST

START_TEST(test_mu_log) {
  for (double x = 0.01; x < 2.0; x += 0.01) {
    ck_assert_ldouble_eq_tol(mu_log(x), log(x), MU_EPS6);
  }

  for (size_t i = 0; i < sizeof(test_cases) / sizeof(test_cases[0]); ++i) {
    if (test_cases[i] > MU_EPS20)
      ck_assert_ldouble_eq_tol(mu_log(test_cases[i]), log(test_cases[i]),
                               MU_EPS6);
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(MU_EPS20, MU_E10);
    ck_assert_ldouble_eq_tol(mu_log(x), log(x), MU_EPS6);
  }

  for (int i = 0; i < 1000; ++i) {
    double x = rand_double(-MU_E10, -MU_EPS20);
    ck_assert_ldouble_nan(mu_log(x));
  }

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
  init_random();
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
