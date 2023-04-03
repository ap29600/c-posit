#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "testing_utils.h"
double random_double(int64_t max_exponent, bool include_negatives) {
	return 1.0
	* (1 - 2 * (rand() % (1 + include_negatives)))
	* (1.0 + (double)rand() / (double)RAND_MAX)
	* pow(2.0, (rand() % (2 * max_exponent)) - max_exponent);
}

int64_t random_int(int64_t max_modulus, bool include_negatives) {
	return 0
	+ (rand() % ((1+include_negatives) * max_modulus))
	- include_negatives * max_modulus;
}

double absf(double d) { return d >= 0 ? d : -d; }

// Note: it only makes sense for these tests to use BW <= 56, as above
// this value posit is not uniformly less precise than double, and
// results can't be verified for correctness.
#define POSIT_BW 32
#define POSIT_SHORT pos
#define POSIT_T posit_t
const uint64_t test_size = 10000;
const uint64_t p_int_max = 1LL << (4 * (POSIT_BW - 3) / 5);
const uint64_t p_expt_max =  4 * (POSIT_BW - 4);
#define POSIT_IMPLEMENTATION
#include "posit.h"


bool integer_conversion(void) {
	test t = {0};

	for (uint64_t i = 0; i < test_size; i++) {
		int64_t n = random_int(p_int_max, true);
		posit_t p = pos_from_i64(n);
		double d = pos_to_double(p);
		expect(&t, d == n);
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}

bool double_conversion(void) {
	test t = {0};

	for (uint64_t i = 0; i < test_size; i++) {
		double n = random_double(p_expt_max, true);
		posit_t p = pos_from_double(n);

		double low = pos_to_double(pos_prev(p));
		double rounded = pos_to_double(p);
		double high = pos_to_double(pos_next(p));

		expect(&t, absf(rounded - n) <= absf(high - n));
		expect(&t, absf(rounded - n) <= absf(low - n));
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}

bool posit_product(void) {
	test t = {0};

	for (uint64_t i = 0; i < test_size; ++i) {
		int64_t n = random_int(sqrt(p_int_max), true);
		int64_t m = random_int(sqrt(p_int_max), true);

		posit_t p = pos_from_i64(n);
		posit_t q = pos_from_i64(m);

		double d = pos_to_double(pos_mul(p, q));
		expect(&t, d == n * m);
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}

bool posit_unary_internal(
	const char *func,
	posit_t posit_fn(posit_t),
	double double_fn(double),
	uint64_t max_expt,
	bool include_negatives
) {
	test t = {.func = func};

	for (uint64_t i = 0; i < test_size; ++i) {
		double d = random_double(max_expt, include_negatives);
		posit_t p = pos_from_double(d);
		posit_t r = posit_fn(p);
		double target = double_fn(pos_to_double(p));

		double low     = pos_to_double(pos_prev(r));
		double rounded = pos_to_double(r);
		double high    = pos_to_double(pos_next(r));

		expect(&t, absf(rounded - target) <= absf(high - target));
		expect(&t, absf(rounded - target) <= absf(low  - target));
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}
#define posit_unary(p, d, m, i) posit_unary_internal(#p, (p), (d), (m), (i))

bool posit_binary_internal(
	const char *func,
	posit_t posit_fn(posit_t, posit_t),
	double double_fn(double, double),
	uint64_t max_expt,
	bool include_negatives
) {
	test t = {.func = func};

	for (uint64_t i = 0; i < test_size; ++i) {
		double d = random_double(max_expt, include_negatives);
		double b = random_double(max_expt, include_negatives);

		posit_t p = pos_from_double(d);
		posit_t q = pos_from_double(b);

		posit_t r = posit_fn(p, q);
		double target = double_fn(pos_to_double(p), pos_to_double(q));

		double low     = pos_to_double(pos_prev(r));
		double rounded = pos_to_double(r);
		double high    = pos_to_double(pos_next(r));

		expect(&t, absf(rounded - target) <= absf(high - target));
		expect(&t, absf(rounded - target) <= absf(low  - target));
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}
#define posit_binary(p, d, m, i) posit_binary_internal(#p, (p), (d), (m), (i))

static double double_inverse(double d)       { return 1 / d; }
static double double_add(double d, double b) { return d + b; }
static double double_mul(double d, double b) { return d * b; }

int main (void) {
	srand(time(0));
	bool ok = true;
	ok &= integer_conversion();
	ok &= double_conversion();
	ok &= posit_product();
	ok &= posit_unary(pos_sqrt, sqrt, p_expt_max, false);
	ok &= posit_unary(pos_inverse, double_inverse, p_expt_max, true);
	ok &= posit_binary(pos_add, double_add, p_expt_max - 1, true);
	ok &= posit_binary(pos_mul, double_mul, p_expt_max / 2, true);
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

