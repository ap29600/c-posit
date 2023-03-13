#include <stdio.h>
#include <inttypes.h>

#define POSIT_BW 32
#define POSIT_ES 2
#define POSIT_IMPLEMENTATION
#include "posit.h"

typedef posit32e2_t p32;
int main (void) {
	quire32_t q = {0};

	q32_add_p32e2(&q, p32e2_from_double(-1));
	q32_debug(&q);
	putchar('\n');

	printf("%g\n", p32e2_to_double(p32e2_from_double(1.1)));
	printf("%g\n", p32e2_to_double(p32e2_from_double(-1.1)));
	q32_add_p32e2(&q, p32e2_from_double(1.1));
	q32_debug(&q);
	putchar('\n');
}

