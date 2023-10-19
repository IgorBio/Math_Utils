#include "math_utils.h"

long int mu_abs(int x) { return x > 0 ? x : -x; }

long double mu_fabs(double x) { return x > 0 ? x : -x; }

long double mu_trunc(double x) { return (x >= 0.0) ? mu_floor(x) : mu_ceil(x); }

long double mu_ceil(double x) {
  if (x >= LLONG_MAX || x <= LLONG_MIN || x != x) {
    return (long double)x;
  }
  long long int_part = (long long)x;
  return (x > int_part) ? (long double)(int_part + 1) : (long double)int_part;
}

long double mu_floor(double x) {
  if (x >= LLONG_MAX || x <= LLONG_MIN || x != x) {
    return (long double)x;
  }
  long long int_part = (long long)x;
  return (x < int_part) ? (long double)(int_part - 1) : (long double)int_part;
}

long double mu_fmod(double x, double y) {
  if (mu_fabs(x) == MU_INF) {
    return MU_NAN;
  } else if (mu_fabs(y) == MU_INF) {
    return x;
  }
  return (long double)(x - mu_trunc(x / y) * y);
}

long double mu_sin(double x) {
  if (x != x || mu_fabs(x) == MU_INF) {
    return MU_NAN;
  }

  x = mu_fmod(x, 2.0 * MU_PI);
  long double res = 0.0;
  long double arg = x;
  unsigned int k = 1;
  int sign = 1;

  while (mu_fabs(arg) > MU_EPS20) {
    res += sign * arg;
    k += 2;
    sign = -sign;
    arg *= (x * x) / (k * (k - 1));
  }

  return res;
}

long double mu_cos(double x) {
  if (x != x || mu_fabs(x) == MU_INF) {
    return MU_NAN;
  }

  x = mu_fmod(x, 2.0 * MU_PI);
  long double res = 0.0;
  long double arg = 1.0;
  unsigned int k = 0;
  int sign = 1;

  while (mu_fabs(arg) > MU_EPS20) {
    res += sign * arg;
    k += 2;
    sign = -sign;
    arg *= (x * x) / (k * (k - 1));
  }

  return res;
}

long double mu_tan(double x) { return mu_sin(x) / mu_cos(x); }

long double mu_asin(double x) {
  if (x != x || mu_fabs(x) > 1.0) {
    return MU_NAN;
  } else if (mu_fabs(x) == 1.0) {
    return MU_PI / 2 * x;
  }

  long double res = x;
  long double arg = x;
  unsigned int k = 1;

  while (mu_fabs(arg) > MU_EPS20) {
    arg *= x * x * (2 * k - 1) * (2 * k - 1) / ((2 * k + 1) * 2 * k);
    res += arg;
    k++;
  }

  return res;
}

long double mu_acos(double x) {
  if (x != x || mu_fabs(x) > 1.0) {
    return MU_NAN;
  }

  return MU_PI / 2 - mu_asin(x);
}

long double mu_atan(double x) {
  if (mu_fabs(x) == MU_INF) {
    return (x > 0) ? MU_PI / 2 : -MU_PI / 2;
  }

  return mu_asin(x / mu_sqrt(1.0 + x * x));
}

long double mu_sqrt(double x) {
  if (x < 0 || x != x || mu_fabs(x) == MU_INF) {
    return MU_NAN;
  }

  long double res = 4.0;
  long double arg = 0.0;

  while (mu_fabs(res - arg) > MU_EPS20) {
    arg = res;
    res = 0.5 * (arg + x / arg);
  }

  return res;
}

long double mu_pow(double base, double exp) {
  if (exp == 0.0 || base == 1.0) {
    return 1.0;
  }

  if (base == 0.0) {
    return (exp != exp) ? MU_NAN : (exp > 0) ? 0.0 : MU_INF;
  }

  if (exp == 1.0) {
    return base;
  }

  if (exp == -1.0) {
    return 1.0 / base;
  }

  if (mu_fabs(base) == MU_INF && (long long)exp != exp) {
    return (exp != exp) ? MU_NAN : (exp > 0) ? MU_INF : 0.0;
  }

  if (base < 0 && (long long)exp != exp) {
    if (mu_fabs(exp) == MU_INF) {
      if (base > -1.0) {
        if (exp == -MU_INF) {
          return MU_INF;
        }
        return 0.0;
      }
      if (base == -1.0) {
        return 1.0;
      }
      if (exp == -MU_INF) {
        return 0.0;
      }
      return MU_INF;
    }
    return MU_NAN;
  }

  if (base != base || exp != exp) {
    return MU_NAN;
  }

  if ((long long)exp == exp) {
    long long int_exp = (long long)exp;
    long double res = 1.0;
    int negative = (int_exp < 0) ? 1 : 0;
    int_exp = mu_fabs(int_exp);

    while (int_exp > 0) {
      if (int_exp & 1) {
        res *= base;
      }
      base *= base;
      int_exp >>= 1;
    }

    return (negative) ? 1.0 / res : res;
  }

  return mu_exp(exp * mu_log(base));
}

long double mu_exp(double x) {
  if (x != x) {
    return MU_NAN;
  }
  if (x == MU_INF) {
    return MU_INF;
  }
  if (x == -MU_INF) {
    return 0.0;
  }

  int negative = 0;
  if (x < 0) {
    x = -x;
    negative = 1;
  }

  long double res = 1.0;
  long double arg = 1.0;
  unsigned int k = 1;

  while (mu_fabs(arg) > MU_EPS20) {
    arg *= x / k;
    res += arg;
    k++;
  }

  res = (negative) ? 1.0 / res : res;

  return res;
}

long double mu_log(double x) {
  if (x == 0.0) {
    return -MU_INF;
  }
  if (x < 0 || x != x) {
    return MU_NAN;
  }
  if (x == MU_INF) {
    return MU_INF;
  }

  long double res = 0.0;
  while (x >= 2.0) {
    x /= 2.0;
    res += MU_LN2;
  }

  x -= 1.0;
  long double arg = x;
  unsigned int k = 1;

  while (mu_fabs(arg) > MU_EPS20) {
    res += arg;
    k++;
    arg *= -x * (k - 1) / k;
  }

  return res;
}
