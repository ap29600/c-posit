#include <stdio.h>
#include <math.h>

// the default bit width is 32
// #define POSIT_BW 32
#define POSIT_IMPLEMENTATION
#include "posit.h"

int main (void) {
	quire32_t q = {0};

	// biggest positive value
	posit32_t big = p32_from_double(pow(2.0, 116));
	// smallest positive value
	posit32_t small = p32_from_double(pow(2.0, -116));
	printf("sums into a quire32 perfectly cancel out:\n");
	q32_add_p32(&q, big);
	q32_add_p32(&q, small);
	q32_add_p32(&q, p32_negate(big));
	q32_add_p32(&q, p32_negate(small));

	q32_debug(&q);
	putchar('\n');
	putchar('\n');

	posit32_t pi = p32_from_double(3.141592654);
	printf( "2 * pi fused multiply-add into a quire32:\n");
	q32_fmadd_p32(&q, pi, p32_from_double(2));
	q32_debug(&q);
	putchar('\n');
}
