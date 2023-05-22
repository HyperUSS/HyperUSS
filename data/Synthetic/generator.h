#include <bits/stdc++.h>

using namespace std;

int possion(double lambda);
double exponential(double lambda);
int zipf(double alpha, int n);
double rand_val(int seed);

int possion(double lambda) {
  double x = -1, u;
  double log1, log2;
  log1 = 0;
  log2 = -lambda;
  do {
    u = rand_val(0);
    log1 += log(u);
    ++x;
  } while (log1 >= log2);
  return x;
}

double exponential(double lambda) {
  double x = 0.0;
  while (true) {
    x = (double)rand() / (double)RAND_MAX;
    if (x != 1) {
      break;
    }
  }
  x = (-1.0 / lambda) * log(1 - x);
  return x;
}

int zipf(double alpha, int n, bool init = false) {
  static double c = 0;      // Normalization constant
  static double *sum_probs; // Pre-calculated sum of probabilities
  double z;                 // Uniform random number (0 < z < 1)
  int zipf_value;           // Computed exponential value to be returned
  int i;                    // Loop counter
  int low, high, mid;       // Binary-search bounds

  // Compute normalization constant on first call only
  if (init && c == 0) {
    for (i = 1; i <= n; i++) {
      c = c + (1.0 / pow((double)i, alpha));
    }
    c = 1.0 / c;

    sum_probs = (double *)malloc((n + 1) * sizeof(*sum_probs));
    sum_probs[0] = 0;
    for (i = 1; i <= n; i++) {
      sum_probs[i] = sum_probs[i - 1] + c / pow((double)i, alpha);
    }
  }

  // Pull a uniform random number (0 < z < 1)
  do {
    z = rand_val(0);
  } while ((z == 0) || (z == 1));

  // Map z to the value
  low = 1, high = n, mid;
  while (low <= high) {
    mid = floor((low + high) / 2);
    if (sum_probs[mid] >= z && sum_probs[mid - 1] < z) {
      zipf_value = mid;
      break;
    } else if (sum_probs[mid] >= z) {
      high = mid - 1;
    } else {
      low = mid + 1;
    }
  }

  if (!((zipf_value >= 1) && (zipf_value <= n))) {
    std::cout << zipf_value << " " << low << " " << high << std::endl;
  }

  // Assert that zipf_value is between 1 and N
  assert((zipf_value >= 1) && (zipf_value <= n));

  return (zipf_value);
}

double rand_val(int seed) {
  const long a = 16807;      // Multiplier
  const long m = 2147483647; // Modulus
  const long q = 127773;     // m div a
  const long r = 2836;       // m mod a
  static long x = 7777;      // Random int value
  long x_div_q;              // x divided by q
  long x_mod_q;              // x modulo q
  long x_new;                // New x value

  // Set the seed if argument is non-zero and then return zero
  if (seed > 0) {
    x = seed;
    return (0.0);
  }

  // RNG using integer arithmetic
  x_div_q = x / q;
  x_mod_q = x % q;
  x_new = (a * x_mod_q) - (r * x_div_q);
  if (x_new > 0) {
    x = x_new;
  } else {
    x = x_new + m;
  }
  // Return a random value between 0.0 and 1.0

  return ((double)x / m);
}
