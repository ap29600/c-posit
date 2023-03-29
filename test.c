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

#define POSIT_BW 32
const uint64_t test_size = 100000;
const uint64_t p_int_max = 1LL << (4 * (POSIT_BW - 3) / 5);
const uint64_t p_expt_max =  4 * (POSIT_BW - 4);
#define POSIT_IMPLEMENTATION
#include "posit.h"


bool integer_conversion(void) {
	test t = {0};

	for (uint64_t i = 0; i < test_size; i++) {
		int64_t n = random_int(p_int_max, true);
		posit32_t p = p32_from_i64(n);
		double d = p32_to_double(p);
		expect(&t, d == n);
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}

bool double_conversion(void) {
	test t = {0};

	for (uint64_t i = 0; i < test_size; i++) {
		double n = random_double(p_expt_max, true);
		posit32_t p = p32_from_double(n);

		double low = p32_to_double(p32_prev(p));
		double rounded = p32_to_double(p);
		double high = p32_to_double(p32_next(p));

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

		posit32_t p = p32_from_i64(n);
		posit32_t q = p32_from_i64(m);

		double d = p32_to_double(p32_mul(p, q));
		expect(&t, d == n * m);
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}

bool posit_unary_internal(
	const char *func,
	posit32_t posit_fn(posit32_t),
	double double_fn(double),
	uint64_t max_expt,
	bool include_negatives
) {
	test t = {.func = func};

	for (uint64_t i = 0; i < test_size; ++i) {
		double d      = random_double(max_expt, include_negatives);
		posit32_t p   = p32_from_double(d);
		posit32_t r   = posit_fn(p);
		double target = double_fn(p32_to_double(p));

		double low     = p32_to_double(p32_prev(r));
		double rounded = p32_to_double(r);
		double high    = p32_to_double(p32_next(r));

		expect(&t, absf(rounded - target) <= absf(high - target));
		expect(&t, absf(rounded - target) <= absf(low  - target));
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}
#define posit_unary(p, d, m, i) posit_unary_internal(#p, (p), (d), (m), (i))

bool posit_binary_internal(
	const char *func,
	posit32_t posit_fn(posit32_t, posit32_t),
	double double_fn(double, double),
	uint64_t max_expt,
	bool include_negatives
) {
	test t = {.func = func};

	for (uint64_t i = 0; i < test_size; ++i) {
		double d      = random_double(max_expt, include_negatives);
		double b      = random_double(max_expt, include_negatives);

		posit32_t p   = p32_from_double(d);
		posit32_t q   = p32_from_double(b);

		posit32_t r   = posit_fn(p, q);
		double target = double_fn(p32_to_double(p), p32_to_double(q));

		double low     = p32_to_double(p32_prev(r));
		double rounded = p32_to_double(r);
		double high    = p32_to_double(p32_next(r));

		expect(&t, absf(rounded - target) <= absf(high - target));
		expect(&t, absf(rounded - target) <= absf(low  - target));
	}

	report(&t);
	return t.ran_tests == t.passed_tests;
}
#define posit_binary(p, d, m, i) posit_binary_internal(#p, (p), (d), (m), (i))

static double double_reciprocal (double d) { return 1 / d; }
static double double_add(double d, double b) { return d + b; }
static double double_mul(double d, double b) { return d * b; }

int main (void) {
	srand(time(0));
	bool ok = true;
	ok &= integer_conversion();
	ok &= double_conversion();
	ok &= posit_product();
	ok &= posit_unary(p32_sqrt, sqrt, p_expt_max, false);
	ok &= posit_unary(p32_reciprocal, double_reciprocal, p_expt_max, true);
	ok &= posit_binary(p32_add, double_add, p_expt_max - 1, true);
	ok &= posit_binary(p32_mul, double_mul, p_expt_max / 2, true);
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

