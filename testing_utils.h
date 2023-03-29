#pragma once

#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

typedef struct {
	uint64_t ran_tests;
	uint64_t passed_tests;
	bool last_test;
	const char *func;
} test;

#define expect(t, expr) expect_internal((t), (expr), # expr, __FILE__, __LINE__, __func__)

#define report(t) test_report_internal((t), __FILE__, __func__)
#define len(a) (sizeof(a) / sizeof(*a))

static inline bool expect_internal(
	test *t,
	bool condition,
	const char *message,
	const char *file,
	uint64_t line,
	const char *func
) {
	t->ran_tests += 1;
	t->last_test = condition;
	if (condition) {
		t->passed_tests += 1;
	} else {
		func = t->func ? t->func : func;
		fprintf(
				stderr,
				"\033[38;5;1m%s:%zu: [test %s]: check `%s` failed.\033[39m\033[49m\n",
				file,
				line,
				func,
				message
		);
	}
	return condition;
}

static inline void test_report_internal(test *t, const char *file, const char *func) {
	char c = 0;
	if (t->passed_tests == t->ran_tests) {
		c = '2';
	} else if (t->passed_tests >  t->ran_tests / 2) {
		c = '3';
	} else {
		c = '1';
	}

	func = t->func ? t->func : func;
	fprintf(
			stderr,
			"\033[38;5;%cm[testing %s/%s]: %zu/%zu checks passed.\033[39m\033[49m\n",
			c,
			file,
			func,
			t->passed_tests,
			t->ran_tests
	);
}
