#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <x86intrin.h>

//=========== CONSTANTS ===============//

#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)
#define CONCAT3(a, b, c) CONCAT(a,CONCAT(b, c))
#define CONCAT4(a, b, c, d) CONCAT(a,CONCAT3(b, c, d))
#define CONCAT5(a, b, c, d, e) CONCAT(a,CONCAT4(b, c, d, e))

#define BITRANGE(a, b) ((uint64_t)((__uint128_t)1 << (a)) - (1ULL << (b)))
#define POSIT_UNPACK(bits) CONCAT(POSIT_SHORT,_unpack)(bits)

#define UNPACKED_POSIT UNPACKED_POSIT

#ifndef POSIT_BW
 #define POSIT_BW 32
#endif

#define POSIT_ES 2
#define QUIRE_BW (16*POSIT_BW)

#define POSIT_SHORT CONCAT(p,POSIT_BW)
#define QUIRE_SHORT CONCAT(q,POSIT_BW)

#define POSIT_T CONCAT3(posit,POSIT_BW,_t)
#define QUIRE_T CONCAT3(quire,POSIT_BW,_t)

#ifndef POSIT_BACKING_T
 #if   POSIT_BW <= 8
  #define POSIT_BACKING_T uint8_t
  #define POSIT_EXTENDED_T uint16_t
 #elif POSIT_BW <= 16
  #define POSIT_BACKING_T uint16_t
  #define POSIT_EXTENDED_T uint32_t
 #elif POSIT_BW <= 32
  #define POSIT_BACKING_T uint32_t
  #define POSIT_EXTENDED_T uint64_t
 #elif POSIT_BW <= 64
  #define POSIT_BACKING_T uint64_t
  #define POSIT_EXTENDED_T __uint128_t
 #else
  #error "posit size is too large: the maximum size is 64."
 #endif
#endif

#ifndef QUIRE_T
 #if   POSIT_BW == 8
  #define QUIRE_T quire8_t
 #elif POSIT_BW == 16
  #define QUIRE_T quire16_t
 #elif POSIT_BW == 32
  #define QUIRE_T quire32_t
 #elif POSIT_BW == 64
  #define QUIRE_T quire64_t
 #else
  #define QUIRE_T CONCAT3(quire_,POSIT_BW,_t)
 #endif
#endif

#if POSIT_BW == 64
 #define POSIT_MASK (-1ULL)
#else
 #define POSIT_MASK ((1ULL << POSIT_BW) - 1)
#endif

#define POSIT_NAR (1ULL << (POSIT_BW - 1))
#define QUIRE_WORDS ((QUIRE_BW - 1) / 64 + 1)
#define QUIRE_SIGN_MASK (1ULL << ((QUIRE_BW - 1) % 64))

//=========== CONSTANTS FOR IEEE754 =============//

#define DOUBLE_EXPT_BIAS     (-1023LL)
#define DOUBLE_BW 64
#define DOUBLE_ES 11
#define DOUBLE_MANTISSA_BW   52
#define DOUBLE_MANTISSA_MASK ((1ULL << DOUBLE_MANTISSA_BW) - 1)
#define DOUBLE_EXPT_MASK     ((1ULL << (DOUBLE_BW - 1)) - (1ULL << DOUBLE_MANTISSA_BW))

//============ TYPE DEFINITIONS ===============//

typedef struct {POSIT_BACKING_T bits;} POSIT_T;
typedef struct {
	uint64_t words[QUIRE_WORDS];
} QUIRE_T __attribute__((aligned (128)));

//========== FUNCTION DECLARATIONS ============//

POSIT_T CONCAT(POSIT_SHORT,_from_i64)    (int64_t i);
POSIT_T CONCAT(POSIT_SHORT,_from_double) (double  d);
double  CONCAT(POSIT_SHORT,_to_double)   (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_negate)      (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_reciprocal)  (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_next)        (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_prev)        (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_sqrt)        (POSIT_T p);

POSIT_T CONCAT(POSIT_SHORT,_mul) (POSIT_T p1, POSIT_T p2);
POSIT_T CONCAT(POSIT_SHORT,_add) (POSIT_T p1, POSIT_T p2);

void CONCAT3(QUIRE_SHORT,_add_,POSIT_SHORT)  (QUIRE_T *q, POSIT_T p);
void CONCAT3(QUIRE_SHORT,_fmadd_,POSIT_SHORT)(QUIRE_T *q, POSIT_T p1, POSIT_T p2);
void CONCAT(QUIRE_SHORT,_debug)              (QUIRE_T const *q);
int  CONCAT(QUIRE_SHORT,_is_nar)             (QUIRE_T const *q);
int  CONCAT(QUIRE_SHORT,_is_zero)            (QUIRE_T const *q);
void CONCAT(QUIRE_SHORT,_set_nar)            (QUIRE_T *q);
void CONCAT(QUIRE_SHORT,_set_zero)           (QUIRE_T *q);

//================= UTILITIES ===================//
// common part shared by different implementation:
// mainly utilities to handle IEEE754 doubles

#ifdef POSIT_IMPLEMENTATION
#ifndef POSIT_IMPLEMENTATION_GUARD
#define POSIT_IMPLEMENTATION_GUARD

struct unpacked_double {
	uint64_t sign;
	int64_t  exponent;
	uint64_t mantissa;
};

static inline __attribute__((always_inline))
int64_t minll(int64_t a, int64_t b) { return a < b ? a : b; }

static inline __attribute__((always_inline))
int64_t maxll(int64_t a, int64_t b) { return a < b ? b : a; }


static inline __attribute__((always_inline))
int64_t clz8(uint8_t a) { return a ? __builtin_clz(a) - 24 : 8 ; }

static inline __attribute__((always_inline))
int64_t clz16(uint16_t a) { return a ? __builtin_clz(a) - 16 : 16; }

static inline __attribute__((always_inline))
int64_t clz32(uint32_t a) { return a ? __builtin_clz(a) : 32; }

static inline __attribute__((always_inline))
int64_t clz64(uint64_t a) { return a ? __builtin_clzll(a) : 64; }

static inline __attribute__((always_inline))
int64_t clz128(__uint128_t a) { return (a >> 64) ? __builtin_clzll(a >> 64) : 64 + clz64(a); }


static inline __attribute__((always_inline))
int64_t absll(int64_t a) { return a < 0 ? -a : a; }

static inline __attribute__((always_inline))
uint64_t extract_bits(uint64_t bits, uint64_t offset, uint64_t size) {
	return (bits << offset >> (64 - size)) & ((1ULL << size) - 1);
}

#define clz(x) _Generic((x), \
	__uint128_t: clz128,   \
	uint64_t:    clz64,    \
	uint32_t:    clz32,    \
	uint16_t:    clz16,    \
	uint8_t:     clz8      \
)(x)

static inline __attribute__((always_inline))
struct unpacked_double unpack_IEEE754_double(double number) {
	union {uint64_t bits; double number;} alias;
	alias.number = number;
	struct unpacked_double result = {
		.sign     = extract_bits(alias.bits, 0, 1),
		.exponent = extract_bits(alias.bits, 1, 11) + DOUBLE_EXPT_BIAS,
		.mantissa = extract_bits(alias.bits, 12, 52),
	};
	return result;
}

#endif

//========= FUNCTION IMPLEMENTATIONS ==========//

static inline __attribute__((always_inline))
POSIT_T CONCAT(POSIT_SHORT,_build)(uint64_t sign, int64_t exponent, POSIT_BACKING_T fraction) {

	if (sign) { exponent = - exponent - 1; }

	uint64_t e = exponent & BITRANGE(POSIT_ES, 0);
	int64_t  r = exponent >> POSIT_ES;
	uint64_t R = r >= 0;
	uint64_t k = minll(POSIT_BW - 2, absll(r) + (r >= 0));

	POSIT_BACKING_T b = 0
		| (POSIT_BACKING_T)!R << (POSIT_BW - 2)
		| e << (POSIT_BW - 2 - POSIT_ES)
		| fraction;

	POSIT_BACKING_T value = 0
		| (POSIT_BACKING_T)sign << (POSIT_BW - 1)
		| R * BITRANGE(POSIT_BW - 1, POSIT_BW - 1 - k)
		| b >> k;

	uint64_t roundup = (b >> (k - 1)) & 1 & (k < POSIT_BW - 2);
	return (POSIT_T){value + roundup};
}

POSIT_T CONCAT(POSIT_SHORT,_from_i64) (int64_t i) {
	uint64_t u = absll(i);
	uint64_t lz = clz(u);
	
	if(lz == 64) { return (POSIT_T){0}; }

	int64_t exponent = 64 - lz - 1;
	uint64_t fraction = u << lz >> (64 - POSIT_BW + POSIT_ES + 1);
	fraction &= BITRANGE(POSIT_BW - POSIT_ES - 2, 0);

	POSIT_T result = CONCAT(POSIT_SHORT,_build)(0, exponent, fraction);
	return i >= 0 ? result : CONCAT(POSIT_SHORT,_negate)(result);
}

// TODO: this doesn't work for subnormal values: we should adjust the exponent
// and shift the mantissa for those; however they are outside the representable
// range for most ES values and will be rounded to minPos anyway.
POSIT_T CONCAT(POSIT_SHORT,_from_double)(double val) {
	if (val == 0) { return (POSIT_T){0}; }
	if (!isfinite(val)) { return (POSIT_T){POSIT_NAR}; }
	struct unpacked_double parts = unpack_IEEE754_double(val);
	uint64_t fraction = parts.mantissa;
	int64_t exponent = parts.exponent;

	// 2's complement and adjust exponent
	if (parts.sign) {
		exponent += !fraction;
		fraction = -fraction & BITRANGE(52, 0);
	}
	fraction <<= 12;
	fraction >>= 64 - POSIT_BW + POSIT_ES + 2;

	return CONCAT(POSIT_SHORT,_build)(parts.sign, exponent, fraction);
}

struct UNPACKED_POSIT {
	uint8_t s;
	int8_t  r;
	uint8_t k;
	uint8_t r0;
	uint8_t e;
	POSIT_BACKING_T f;
};

static inline __attribute__((always_inline))
struct UNPACKED_POSIT POSIT_UNPACK(POSIT_T p) {
	struct UNPACKED_POSIT result = {0};
	uint64_t s  = p.bits >> (POSIT_BW - 1);
	uint64_t r0 = (p.bits >> (POSIT_BW - 2)) & 1;
	uint64_t k  = __builtin_clzll(((((uint64_t)p.bits ^ -r0) << 1) & POSIT_MASK) | 1) - (64 - POSIT_BW);
	uint64_t r  = r0 ? k - 1 : -k;
	uint64_t e  = ((uint64_t)p.bits << k >> (POSIT_BW - 2 - POSIT_ES)) & BITRANGE(POSIT_ES, 0);
	uint64_t f  = ((uint64_t)p.bits << k << (POSIT_ES + 2)) & BITRANGE(POSIT_BW, 0);

	result.s = s;
	result.r0 = r0;
	result.k = k;
	result.r = r;
	result.e = e;
	result.f = f;
	return result;
}


double CONCAT(POSIT_SHORT,_to_double)(POSIT_T val) {
	if (val.bits == 0) { return 0.0; }
	if (val.bits == POSIT_NAR) { return NAN; }

	struct UNPACKED_POSIT parts = POSIT_UNPACK(val);

	uint64_t posit_fraction = 1
		* (1 - (parts.s << 1))
		* ((uint64_t)parts.f << (64 - POSIT_BW) >> 1);
	uint64_t mantissa = posit_fraction >> (64 - 52 - 1);
	uint64_t roundup = (posit_fraction >> (64 - 52 - 2)) & 1;

	int64_t exponent = 1
		* (1 - (parts.s << 1))
		* (parts.r * (1 << POSIT_ES) + parts.e + (parts.s & (parts.f > 0)));
	
	mantissa += roundup;
	mantissa &= DOUBLE_MANTISSA_MASK;
	exponent += roundup * (mantissa == 0);

	union{uint64_t bits; double value;} result;
	result.bits = 0
		| (uint64_t)parts.s << 63
		| (uint64_t)(exponent - DOUBLE_EXPT_BIAS) << DOUBLE_MANTISSA_BW
		| mantissa;
	return result.value;
}

POSIT_T CONCAT(POSIT_SHORT,_negate)(POSIT_T val) {
	val.bits = - val.bits & POSIT_MASK;
	return val;
}

POSIT_T CONCAT(POSIT_SHORT,_floor)(POSIT_T val) {
	// as an unsigned integer < posit(1)
	//   ==> as a posit it's in the range [0, 1)
	if (val.bits < (1ULL << (POSIT_BW - 2))) { return (POSIT_T) {0}; }

	// as an unsigned integer >= posit(-1)
	//   ==> as a posit it's in the range [-1, 0)
	if (val.bits >= (-(1ULL << (POSIT_BW - 2)) & POSIT_MASK)) {
		return (POSIT_T) {-(1ULL << (POSIT_BW - 2)) & POSIT_MASK};
	}

	// this handles NaR as well, because in that case we get k = POSIT_BW - 1
	struct UNPACKED_POSIT parts = POSIT_UNPACK(val);
	uint64_t offset = 0
		+ parts.k * ((1ULL << POSIT_ES) + 1)
		+ (parts.s
			? -(uint64_t)parts.e - 1
			: parts.e - (1ULL << POSIT_ES));

	// mask out the bottom bits
	val.bits &= BITRANGE(
			POSIT_BW,
			maxll(0, POSIT_BW - POSIT_ES - 2 - offset)
	);

	return val;
}

// TODO: write a dedicated function instead of this
POSIT_T CONCAT(POSIT_SHORT,_ceil)(POSIT_T val) {
	POSIT_T n = CONCAT(POSIT_SHORT,_negate)(val);
	POSIT_T cn = CONCAT(POSIT_SHORT,_floor)(n);
	return CONCAT(POSIT_SHORT,_negate)(cn);
}

POSIT_T CONCAT(POSIT_SHORT,_next)(POSIT_T val) {
	return (POSIT_T){.bits = (val.bits + 1) & POSIT_MASK};
}

POSIT_T CONCAT(POSIT_SHORT,_prev)(POSIT_T val) {
	return (POSIT_T){.bits = (val.bits - 1) & POSIT_MASK};
}

POSIT_T CONCAT(POSIT_SHORT,_reciprocal)(POSIT_T p) {
	if (p.bits == POSIT_NAR) { return p; }
	if (p.bits == 0) { return (POSIT_T){POSIT_NAR}; }

	struct UNPACKED_POSIT parts = POSIT_UNPACK(p);

	// TODO: figure out the correct logic for negative numbers.
	if (parts.s) {
		POSIT_T n = CONCAT(POSIT_SHORT,_negate)(p);
		POSIT_T rn = CONCAT(POSIT_SHORT,_reciprocal)(n);
		return CONCAT(POSIT_SHORT,_negate)(rn);
	}

	int64_t exponent = 1
		* (1 - (parts.s << 1))
		* (parts.r * (1 << POSIT_ES) + parts.e)
		* -1;

	POSIT_EXTENDED_T one = (POSIT_EXTENDED_T)1ULL << (2 * POSIT_BW - 1);
	POSIT_BACKING_T xf = POSIT_NAR | parts.f >> 1;
	POSIT_BACKING_T quot = (one / (POSIT_EXTENDED_T)xf) & POSIT_MASK;

	exponent -= quot > 0;
	quot >>= 2 + POSIT_ES - 1;
	quot &= BITRANGE(POSIT_BW - POSIT_ES - 2, 0);

	return CONCAT(POSIT_SHORT,_build)(parts.s, exponent, quot);
}

POSIT_T CONCAT(POSIT_SHORT,_sqrt) (POSIT_T p) {
	if (p.bits == 0) { return p; }
	// negative numbers and NaR
	if (p.bits >= POSIT_NAR) { return (POSIT_T){POSIT_NAR}; }

	struct UNPACKED_POSIT parts = POSIT_UNPACK(p);

	int64_t old_exponent = parts.r * (1 << POSIT_ES) + parts.e;
	int64_t new_exponent = (old_exponent + 1) >> 1;
	int64_t exponent_rem = old_exponent & 1;


	POSIT_BACKING_T f = 0
		| (POSIT_BACKING_T)1 << (POSIT_BW - 1)
		| parts.f >> 1;

	POSIT_EXTENDED_T t = (POSIT_EXTENDED_T)f << (POSIT_BW - 1 - exponent_rem);
	POSIT_BACKING_T residual = 0;

	// Newton Iteration
	do {
		f -= residual;
		// TODO: is it better to just do more iterations with the secant method?
		residual = ((POSIT_EXTENDED_T)f * (POSIT_EXTENDED_T)f - t) / f >> 1;
	} while(residual);

	f -= ((POSIT_EXTENDED_T)f * (POSIT_EXTENDED_T)f) > t;

	uint64_t lz = clz(f) - (8 * sizeof(f) - POSIT_BW);
	new_exponent -= lz;
	f <<= lz + 1;
	f >>= POSIT_ES + 2;
	f &= BITRANGE(POSIT_BW - POSIT_ES - 2, 0);

	return CONCAT(POSIT_SHORT,_build)(0, new_exponent, f);
}


POSIT_T CONCAT(POSIT_SHORT,_add) (POSIT_T p1, POSIT_T p2) {
	if (p1.bits == POSIT_NAR || p2.bits == POSIT_NAR) { return (POSIT_T){POSIT_NAR}; }
	if (p1.bits == 0) { return p2; }
	if (p2.bits == 0) { return p1; }

	struct UNPACKED_POSIT parts1 = POSIT_UNPACK(p1);
	struct UNPACKED_POSIT parts2 = POSIT_UNPACK(p2);

	int64_t shift1 = 1
		* (1ULL - (parts1.s << 1))
		* ((1ULL << POSIT_ES) * parts1.r + parts1.e);

	int64_t shift2 = 1
		* (1ULL - (parts2.s << 1))
		* ((1ULL << POSIT_ES) * parts2.r + parts2.e);

	int64_t shift = maxll(shift1, shift2);

	POSIT_EXTENDED_T f1 = 0
		| -(POSIT_EXTENDED_T)parts1.s << (POSIT_BW - shift + shift1)
		| 1ULL << (POSIT_BW - 1 - shift + shift1)
		| (POSIT_EXTENDED_T)parts1.f >> (1 + shift - shift1 + parts1.s);

	POSIT_EXTENDED_T f2 = 0
		| -(POSIT_EXTENDED_T)parts2.s << (POSIT_BW - shift + shift2)
		| 1ULL << (POSIT_BW - 1 - shift + shift2)
		| (POSIT_EXTENDED_T)parts2.f >> (1 + shift - shift2 + parts2.s);

	f1 = (shift - shift1 < POSIT_BW) * f1;
	f2 = (shift - shift2 < POSIT_BW) * f2;

	POSIT_EXTENDED_T f = f1 + f2;
	uint64_t sign = (f >> (2 * POSIT_BW - 1)) & 1;
	uint64_t lz = clz(f ^ -(POSIT_EXTENDED_T)sign) - (sizeof(f) * 8 - 2 * POSIT_BW);

	if (f == 0) { return (POSIT_T){0}; }
	f <<= lz + 1;
	f >>= POSIT_BW + POSIT_ES + 2;
	f &= BITRANGE(POSIT_BW - POSIT_ES - 2, 0);

	int64_t exponent = shift - lz + POSIT_BW;
	return CONCAT(POSIT_SHORT,_build)(sign, exponent, f);
}

POSIT_T CONCAT(POSIT_SHORT,_mul) (POSIT_T p1, POSIT_T p2) {
	if (p1.bits == POSIT_NAR || p2.bits == POSIT_NAR) { return (POSIT_T){POSIT_NAR}; }
	if (p1.bits == 0ULL || p2.bits == 0ULL) { return (POSIT_T) {0}; }

	struct UNPACKED_POSIT parts1 = POSIT_UNPACK(p1);
	struct UNPACKED_POSIT parts2 = POSIT_UNPACK(p2);

	int64_t shift1 = 1
		* (1ULL - (parts1.s << 1))
		* ((1ULL << POSIT_ES) * parts1.r + parts1.e);

	int64_t shift2 = 1
		* (1ULL - (parts2.s << 1))
		* ((1ULL << POSIT_ES) * parts2.r + parts2.e);

	int64_t shift = shift1 + shift2;

	POSIT_EXTENDED_T f1 = 0
		| -(POSIT_EXTENDED_T)parts1.s << POSIT_BW
		| 1ULL << (POSIT_BW - 1)
		| (POSIT_BACKING_T)parts1.f >> (1 + parts1.s);

	POSIT_EXTENDED_T f2 = 0
		| -(POSIT_EXTENDED_T)parts2.s << POSIT_BW
		| 1ULL << (POSIT_BW - 1)
		| (POSIT_BACKING_T)parts2.f >> (1 + parts2.s);

	POSIT_EXTENDED_T f = f1 * f2;
	uint64_t sign = parts1.s ^ parts2.s;
	uint64_t lz = clz(f ^ -(POSIT_EXTENDED_T)sign) - (8 * sizeof(f) - 2 * POSIT_BW);
	f <<= lz + 1;
	f >>= POSIT_BW + POSIT_ES + 2;
	f &= BITRANGE(POSIT_BW - POSIT_ES - 2, 0);
	
	int64_t  exponent = shift - lz + 1;
	return CONCAT(POSIT_SHORT,_build)(sign, exponent, f);
}

size_t CONCAT(POSIT_SHORT,_to_bits)(char *buf, size_t len, POSIT_T p) {
	size_t written = 0;
	size_t i = 0;
	struct UNPACKED_POSIT parts = POSIT_UNPACK(p);

	if (written + 1 >= len || ++i > POSIT_BW) goto ret;
	buf[written++] = '0' + parts.s;

	if (written + 1 >= len) goto ret;
	buf[written++] = ' ';

	for (int j = 0; j < parts.k; ++j) {
		if (written + 1 >= len || ++i > POSIT_BW) goto ret;
		buf[written++] = '0' + parts.r0;
	}
	if (written + 1 >= len || ++i > POSIT_BW) goto ret;
	buf[written++] = '0' + !parts.r0;

	if (written + 1 >= len) goto ret;
	buf[written++] = ' ';
	
	for (int j = 0; j < POSIT_ES; ++j) {
		if (written + 1 >= len || ++i > POSIT_BW) goto ret;
		buf[written++] = '0' + ((parts.e >> (POSIT_ES - j - 1)) & 1);
	}

	if (written + 1 >= len) goto ret;
	buf[written++] = ' ';

	for (int j = 0; j < POSIT_BW - parts.k - 2 - POSIT_ES; ++j) {
		if (written + 1 >= len || ++i > POSIT_BW) goto ret;
		buf[written++] = '0' + ((parts.f >> (POSIT_BW - j - 1)) & 1);
	}

	ret:
	if (written >= len) return written;
	buf[written++] = '\0';

	return written;
}

void CONCAT3(QUIRE_SHORT,_add_,POSIT_SHORT)(QUIRE_T *q, POSIT_T p) {
	if (p.bits == 0ULL) { return; }

	if (p.bits == POSIT_NAR) {
		CONCAT(QUIRE_SHORT,_set_nar)(q);
		return;
	}

	if (CONCAT(QUIRE_SHORT,_is_nar)(q)) { return; }

	struct UNPACKED_POSIT parts = POSIT_UNPACK(p);

	int64_t shift = 0
		+ (8 * POSIT_BW - 16) - 63
		+ (1ULL - (parts.s << 1)) * ((1ULL << POSIT_ES) * parts.r + parts.e);

	uint64_t bits = 0
		| 1ULL << 63
		| (uint64_t)parts.f << (64 - POSIT_BW) >> (1 + parts.s);

	uint64_t fill_word = parts.s * -1ULL;
	
	uint64_t hi_bits = fill_word << (shift & 63);
	if (shift & 63) { hi_bits |= bits >> (64 - (shift & 63)); }
	uint64_t lo_bits = bits << (shift & 63);

	uint64_t carry = 0;

	if (shift >> 6 >= 0) {
		carry = _addcarry_u64(
				carry,
				q->words[shift >> 6],
				lo_bits,
				(unsigned long long*)&q->words[shift >> 6]
		);
	}
	carry = _addcarry_u64(
			carry,
			q->words[(shift >> 6) + 1],
			hi_bits,
			(unsigned long long*)&q->words[(shift >> 6) + 2]
	);
	for (int64_t i = (shift >> 6) + 2; i < QUIRE_WORDS; i++ ) {
		carry = _addcarry_u64(
				carry,
				q->words[i],
				fill_word,
				(unsigned long long*)&q->words[i]
		);
	}

	q->words[QUIRE_WORDS - 1] &= (QUIRE_SIGN_MASK << 1) - 1;
}

void CONCAT3(QUIRE_SHORT,_fmadd_,POSIT_SHORT)(QUIRE_T *q, POSIT_T p1, POSIT_T p2) {
	
	if (p1.bits == 0ULL) { return; }
	if (p2.bits == 0ULL) { return; }

	if (p1.bits == POSIT_NAR || p2.bits == POSIT_NAR) {
		CONCAT(QUIRE_SHORT,_set_nar)(q);
		return;
	}

	if (CONCAT(QUIRE_SHORT,_is_nar)(q)) { return; }

	struct UNPACKED_POSIT parts1 = POSIT_UNPACK(p1);
	struct UNPACKED_POSIT parts2 = POSIT_UNPACK(p2);

	int64_t shift1 = 1
		* (1ULL - (parts1.s << 1))
		* ((1ULL << POSIT_ES) * parts1.r + parts1.e);

	int64_t shift2 = 1
		* (1ULL - (parts2.s << 1))
		* ((1ULL << POSIT_ES) * parts2.r + parts2.e);

	uint64_t fill_word = -((uint64_t)parts1.s ^ (uint64_t)parts2.s);
	int64_t shift = shift1 + shift2 + (8 * POSIT_BW - 16) - 126;

	__uint128_t f1 = 0
		| ((__uint128_t)-(uint64_t) parts1.s) << 64
		| 1ULL << 63
		| (uint64_t)parts1.f << (64 - POSIT_BW) >> (1 + parts1.s);

	__uint128_t f2 = 0
		| ((__uint128_t)-(uint64_t) parts2.s) << 64
		| 1ULL << 63
		| (uint64_t)parts2.f << (64 - POSIT_BW) >> (1 + parts2.s);

	__uint128_t f = f1 * f2;

	uint64_t hi_bits = fill_word << (shift & 63);
	if (shift & 63) { hi_bits |= f >> (128 - (shift & 63)); }
	uint64_t md_bits = f >> (64 - (shift & 63));
	uint64_t lo_bits = f << (shift & 63);

	uint64_t carry = 0;
	if (shift >> 6 >= 0) {
		carry = _addcarry_u64(
				carry,
				q->words[shift >> 6],
				lo_bits,
				(unsigned long long*)&q->words[shift >> 6]
		);
	}
	carry = _addcarry_u64(
			carry,
			q->words[(shift >> 6) + 1],
			md_bits,
			(unsigned long long*)&q->words[(shift >> 6) + 1]
	);
	carry = _addcarry_u64(
			carry,
			q->words[(shift >> 6) + 2],
			hi_bits,
			(unsigned long long*)&q->words[(shift >> 6) + 2]
	);
	for (int64_t i = (shift >> 6) + 3; i < QUIRE_WORDS; i++ ) {
		carry = _addcarry_u64(
				carry,
				q->words[i],
				fill_word,
				(unsigned long long*)&q->words[i]
		);
	}

	q->words[QUIRE_WORDS - 1] &= (QUIRE_SIGN_MASK << 1) - 1;
}

void CONCAT(QUIRE_SHORT,_negate)(QUIRE_T *q) {
	uint64_t carry = 1;
	for(int64_t i = 0; i < QUIRE_WORDS; i++) {
		carry = _addcarry_u64(
				carry,
				~q->words[i],
				0,
				(unsigned long long *)&q->words[i]
		);
	}
	// discard any overflows
	q->words[QUIRE_WORDS - 1] &= (QUIRE_SIGN_MASK << 1) - 1;
}

POSIT_T CONCAT3(QUIRE_SHORT,_to_,POSIT_SHORT)(QUIRE_T const *q) {
	if (q->words[QUIRE_WORDS - 1] & QUIRE_SIGN_MASK) {
		assert(0 && "unimplemented");
	} else {
		int64_t lz = ((QUIRE_BW) % 64)
				? (- 64 + ((QUIRE_BW) % 64))
				: 0;

		for(int i = QUIRE_WORDS - 1; i >= 0; --i) {
			uint64_t lzw = clz(q->words[i]);
			lz += lzw;
			if (lzw != 64) break;
		}
		if (lz == QUIRE_BW) { return (POSIT_T){0}; }

		int64_t offset = QUIRE_BW - lz - 1;
		uint64_t bits = q->words[offset >> 6] << (63 - (offset & 63));

		if ((offset >> 6) > 0 && ((offset & 63) < 63)) {
			bits |= q->words[(offset >> 6) - 1] >> ((offset & 63) + 1);
		}
		bits >>= 64 - POSIT_BW + POSIT_ES + 1;
		bits &= BITRANGE(POSIT_BW - POSIT_ES - 2, 0);

		int64_t exponent = offset - (8 * POSIT_BW - 16);
		return CONCAT(POSIT_SHORT,_build)(0, exponent, bits);
	}
}

void CONCAT(QUIRE_SHORT,_debug) (QUIRE_T const*q) {
	for(uint64_t i = 0; i < 2 * POSIT_BW; ++i) {
		if (i == 4) { putchar('|'); };
		if (i == 4 + POSIT_BW - 2)  { printf("\n        ."); }
		uint64_t j = 2 * POSIT_BW - i - 1;
		uint8_t c = q->words[j / 8] >> (j % 8) * 8;
		printf("%02x", c);
	}
 }

int CONCAT(QUIRE_SHORT,_is_nar)(QUIRE_T const *q) {
	int result = 1;
	for (int i = 0; i < QUIRE_WORDS - 1; ++i) {
		result &= q->words[i] == 0;
	}
	result &= q->words[QUIRE_WORDS - 1] == QUIRE_SIGN_MASK;
	return result;
}

int CONCAT(QUIRE_SHORT,_is_zero)(QUIRE_T const *q) {
	int result = 1;
	for (int i = 0; i < QUIRE_WORDS; ++i) {
		result &= q->words[i] == 0;
	}
	return result;
}

void CONCAT(QUIRE_SHORT,_set_nar)(QUIRE_T *q) {
	for (int i = 0; i < QUIRE_WORDS - 1; ++i) {
		q->words[i] = 0;
	}
	q->words[QUIRE_WORDS - 1] = QUIRE_SIGN_MASK;
}

void CONCAT(QUIRE_SHORT,_set_zero)(QUIRE_T *q) {
	for (int i = 0; i < QUIRE_WORDS - 1; ++i) {
		q->words[i] = 0;
	}
}

//=========== CONSTANT CLEANUP ===============//

#endif

#undef POSIT_UNPACK
#undef UNPACKED_POSIT

#undef POSIT_SHORT
#undef POSIT_T
#undef POSIT_BW
#undef POSIT_ES
#undef POSIT_BACKING_T
#undef POSIT_EXTENDED_T
#undef POSIT_NAR
#undef POSIT_MASK
#undef QUIRE_T

#undef DOUBLE_EXPT_BIAS
#undef DOUBLE_BW
#undef DOUBLE_ES
#undef DOUBLE_MANTISSA_BW
#undef DOUBLE_MANTISSA_MASK
#undef DOUBLE_EXPT_MASK

#undef CONCAT_
#undef CONCAT
#undef CONCAT3
#undef CONCAT4
#undef CONCAT5
