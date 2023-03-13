#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

//=========== CONSTANTS ===============//

#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)
#define CONCAT3(a, b, c) CONCAT(a,CONCAT(b, c))
#define CONCAT4(a, b, c, d) CONCAT(a,CONCAT3(b, c, d))
#define CONCAT5(a, b, c, d, e) CONCAT(a,CONCAT4(b, c, d, e))

#define POSIT_BITS(bits, start, size) CONCAT(POSIT_SHORT,_extract_bits)(bits, start, size)
#define POSIT_UNPACK(bits) CONCAT(POSIT_SHORT,_unpack)(bits)
#define UNPACKED_POSIT UNPACKED_POSIT

#ifndef POSIT_BW
 #define POSIT_BW 32
#endif

#define POSIT_ES 2

#define POSIT_SHORT CONCAT4(p,POSIT_BW,e,POSIT_ES)
#define QUIRE_SHORT CONCAT(q,POSIT_BW)

#define POSIT_T CONCAT5(posit,POSIT_BW,e,POSIT_ES,_t)
#define QUIRE_T CONCAT3(quire,POSIT_BW,_t)

#ifndef POSIT_BACKING_T
 #if   POSIT_BW <= 8
  #define POSIT_BACKING_T uint8_t
 #elif POSIT_BW <= 16
  #define POSIT_BACKING_T uint16_t
 #elif POSIT_BW <= 32
  #define POSIT_BACKING_T uint32_t
 #elif POSIT_BW <= 64
  #define POSIT_BACKING_T uint64_t
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
  #define QUIRE_T quire_t
 #endif
#endif

#if POSIT_BW == 64
 #define POSIT_MASK (-1ULL)
#else
 #define POSIT_MASK ((1ULL << POSIT_BW) - 1)
#endif

#define POSIT_NAR (1ULL << (POSIT_BW - 1))
#define QUIRE_WORDS ((POSIT_BW * 16 - 1) / 64 + 1)
#define QUIRE_SIGN_MASK (1ULL << ((16 * POSIT_BW - 1) % 64))

//=========== CONSTANTS FOR IEEE754 =============//

#define DOUBLE_EXPT_BIAS     (-1023LL)
#define DOUBLE_BW 64
#define DOUBLE_ES 11
#define DOUBLE_MANTISSA_BW   52
#define DOUBLE_MANTISSA_MASK ((1ULL << DOUBLE_MANTISSA_BW) - 1)
#define DOUBLE_EXPT_MASK     ((1ULL << (DOUBLE_BW - 1)) - (1ULL << DOUBLE_MANTISSA_BW))

//============ TYPE DEFINITIONS ===============//

typedef struct {POSIT_BACKING_T bits;} POSIT_T;
typedef struct {uint64_t words[QUIRE_WORDS];} QUIRE_T;

//========== FUNCTION DECLARATIONS ============//

POSIT_T CONCAT(POSIT_SHORT,_from_double)   (double d);
double  CONCAT(POSIT_SHORT,_to_double)     (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_negate)        (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_next)          (POSIT_T p);
POSIT_T CONCAT(POSIT_SHORT,_prev)          (POSIT_T p);
void CONCAT3(QUIRE_SHORT,_add_,POSIT_SHORT)(QUIRE_T *q, POSIT_T p);
void CONCAT(QUIRE_SHORT,_debug)            (QUIRE_T *q);

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
int64_t absll(int64_t a) { return a < 0 ? -a : a; }

static inline __attribute__((always_inline))
uint64_t extract_bits(uint64_t bits, uint64_t offset, uint64_t size) {
	return (bits << offset >> (64 - size)) & ((1ULL << size) - 1);
}

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

// TODO: this doesn't work for subnormal values: we should adjust the exponent
// and shift the mantissa for those; however they are outside the representable
// range for most ES values and will be rounded to minPos anyway.
POSIT_T CONCAT(POSIT_SHORT,_from_double)(double val) {
	POSIT_BACKING_T value = 0;
	if (val == 0) { goto ret; }
	if (!isfinite(val)) { value = POSIT_NAR; goto ret; }
	struct unpacked_double parts = unpack_IEEE754_double(val);
	POSIT_BACKING_T fraction = extract_bits(parts.mantissa, 12, POSIT_BW - 3 - POSIT_ES);

	// 2's complement and adjust exponent
	if (parts.sign) {
		fraction = ((1ULL << (POSIT_BW - 3 - POSIT_ES)) - 1) & -fraction;
		parts.exponent = - parts.exponent - !!fraction;
	}

	uint64_t e = parts.exponent & ((1ULL << POSIT_ES) - 1);
	int64_t  r = parts.exponent >> POSIT_ES;
	uint64_t R = r >= 0;
	uint64_t k = minll(POSIT_BW - 3, absll(r) - (r < 0));

	// TODO: this always rounds down, we should try to round to
	// nearest.
	POSIT_BACKING_T b =
		(POSIT_BACKING_T)!R << (POSIT_BW - 3)
		| e << (POSIT_BW - 3 - POSIT_ES)
		| fraction;

	value =
		parts.sign << (POSIT_BW - 1)
		| R * ((1ULL << (POSIT_BW - 1)) - (1ULL << (POSIT_BW - 2 - k)))
		| b >> k;

	ret:
	return (POSIT_T){value};
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
	uint64_t e  = (p.bits << k >> (POSIT_BW - 2 - POSIT_ES)) & ((1ULL << POSIT_ES) - 1);
	uint64_t f  = (p.bits << k << (POSIT_ES + 2)) & POSIT_MASK;

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

	int64_t exponent;
	uint64_t mantissa;
	if (parts.s) {
		mantissa = - extract_bits((uint64_t)parts.f, 64 - POSIT_BW, 52) & DOUBLE_MANTISSA_MASK;
		exponent = - (int)(1 << POSIT_ES) * parts.r + parts.e - (parts.f > 0);
	} else {
		exponent = (int)(1 << POSIT_ES) * parts.r + parts.e;
		mantissa = extract_bits((uint64_t)parts.f, 64 - POSIT_BW, 52);
	}

	union{uint64_t bits; double value;} result;
	result.bits =
		(uint64_t)parts.s << 63
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
	if (val.bits < (1ULL << (POSIT_BW - 2))) {
		return (POSIT_T) {0};
	}

	// as an unsigned integer >= posit(-1)
	//   ==> as a posit it's in the range [-1, 0)
	if (val.bits >= (-(1ULL << (POSIT_BW - 2)) & POSIT_MASK)) {
		return (POSIT_T) {-(1ULL << (POSIT_BW - 2)) & POSIT_MASK};
	}

	// this handles NaR as well, because in that case we get k = POSIT_BW - 1
	struct UNPACKED_POSIT parts = POSIT_UNPACK(val);
	uint64_t offset =
		parts.k * ((1ULL << POSIT_ES) + 1)
		+ (parts.s
			? -(uint64_t)parts.e - 1
			: parts.e - (1ULL << POSIT_ES));

	POSIT_BACKING_T mask =
		POSIT_MASK
		>> (2 + POSIT_ES)
		>> minll(POSIT_BW - 2 - POSIT_ES, offset);

	val.bits &= ~mask;
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
		q->words[QUIRE_WORDS - 1] = QUIRE_SIGN_MASK;
		for (int i = 0; i < QUIRE_WORDS - 1; ++i) {
			q->words[i] = 0;
		}
		return;
	}

	{
		int is_nar = q->words[QUIRE_WORDS - 1] == QUIRE_SIGN_MASK;
		for (int i = 0; i < QUIRE_WORDS - 1; ++i) {
			is_nar &= q->words[i] == 0;
		}
		if (is_nar) { return; }
	}

	struct UNPACKED_POSIT parts = POSIT_UNPACK(p);

	int64_t shift =
		(8 * POSIT_BW - 16) - 63
		+ (parts.s
			? (-(1LL << POSIT_ES) * parts.r + parts.e - 1)
			: ( (1LL << POSIT_ES) * parts.r + parts.e));

	uint64_t bits =
		(uint64_t)!parts.s << 63
		| (uint64_t)parts.f << (64 - POSIT_BW) >> 1;

	uint64_t fill_word = parts.s ? -1ULL : 0;
	uint64_t carry = 0;
	
	uint64_t lo_bits = bits << (shift & 63);
	uint64_t hi_bits = fill_word;
	if (shift & 63) {
		hi_bits &= ~(-1ULL >> (64 - (shift & 63)));
		hi_bits |= bits >> (64 - (shift & 63));
	}

	if (shift >= 0) {
		uint64_t prev = q->words[shift >> 6];
		uint64_t word = lo_bits;
		uint64_t add = prev + word;
		carry = add < prev;
		q->words[shift >> 6] = add;
	}

	if ((shift & 63) || carry) {
		uint64_t prev = q->words[(shift >> 6) + 1];
		uint64_t word = hi_bits;
		uint64_t add = prev + word + carry;
		q->words[(shift >> 6) + 1] = add;
		carry = fill_word + (add < prev || add < word);
	}

	// propagate the carry
	uint64_t i = (shift >> 6) + 2;
	for(; carry && i < QUIRE_WORDS; ++i) {
		uint64_t word = carry;
		uint64_t prev = q->words[i];
		uint64_t add  = word + prev;
		carry = fill_word + (add < word);
		q->words[i] += word;
	}

	// discard any overflows
	q->words[QUIRE_WORDS - 1] &= (QUIRE_SIGN_MASK << 1) - 1;
}

void CONCAT(QUIRE_SHORT,_debug) (QUIRE_T *q) {
	for(uint64_t i = 0; i < 2 * POSIT_BW; ++i) {
		if (i == 4) { putchar('|'); };
		if (i == 4 + POSIT_BW - 2)  { putchar('.'); }
		uint64_t j = 2 * POSIT_BW - i - 1;
		uint8_t c = q->words[j / 8] >> (j % 8) * 8;
		printf("%02x", c);
	}
}

//=========== CONSTANT CLEANUP ===============//

#endif

#undef POSIT_BITS
#undef POSIT_UNPACK
#undef UNPACKED_POSIT

#undef POSIT_SHORT
#undef POSIT_T
#undef POSIT_BW
#undef POSIT_ES
#undef POSIT_BACKING_T
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
