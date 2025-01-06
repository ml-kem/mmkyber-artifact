//  zp_dist.h
//  === Header: discrete distributions

#ifndef _ZP_DIST_H_
#define _ZP_DIST_H_

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <gmp.h>

//  discrete distribution

typedef struct {
    int64_t x;          //  x at first entry
    size_t eps;         //  target precision, bits
    size_t n;           //  number of entries; range is x .. x+n-1
    size_t v_sz;        //  allocated array size
    mpz_t *v;           //  distribution scaled by 2^eps
} dist_t;

//  intialize operations
int dist_init(size_t eps);

//  Return a new distribution with "n" elements, precision "eps" bits.
dist_t *dist_alloc(size_t n, size_t eps);

//  Make sure distribution "d" has space for "sz" values, set to 0.
dist_t *dist_grow(dist_t *dr, size_t sz, size_t eps);

//  zeroize a vector
dist_t *dist_zero(dist_t *dr, size_t n, size_t eps);

//  set dr to da
dist_t *dist_copy(dist_t *dr, const dist_t *da);

//  make a new copy
dist_t *dist_dup(const dist_t *da);

//  Free distribution "d"
void dist_clear(dist_t *d);

//  print information about distribution "d", labeled "lab"
//  returns approximate variance
double dist_print(const dist_t *d, const char *lab);

//  (debug:) dumps internal vector as a polynomial (suitable for gp/pari)
void vecz_dump(mpz_t *v, size_t n, const char *lab);

//  Scale to sum to 2**eps, prune. If lab != NULL, report size.
size_t dist_sum1(dist_t *d, const char *lab);

//  create a discrete gaussian distribution with given parameters
//  width=false: use std. devation. true: scale sigma by 1/sqrt(2*Pi)
dist_t *dist_gauss( dist_t *dr, const char *mu,
                    const char *sig, bool width, size_t eps);

//  create a uniform distribution in range [a, b]
dist_t *dist_unif(dist_t *dr, int64_t a, int64_t b, size_t eps);

//  centered binomial distribution [-eta, eta]
dist_t *dist_cbd(dist_t *dr, int64_t eta, size_t eps);

//  Psi, the roudning distribution: r - Decompress_q(Compress_d(r))
dist_t *dist_round(dist_t *dr, int64_t d, int64_t q, size_t eps);

//  dbl()-and-round distribution: c_i - [dbl(c_i]_{2^{d_u}}
dist_t *dist_dbl_r(dist_t *dr, int64_t d, int64_t q, size_t eps);

//  negate a distribution: dr = -da
dist_t *dist_neg(dist_t *dr, const dist_t *da);

//  Add (convolute) two distributions "da" and "db", result to "dr"
dist_t *dist_add(dist_t *dr, const dist_t *da, const dist_t *db);

//  Signed multiply of two distributions
dist_t *dist_mul(dist_t *dr, const dist_t *da, const dist_t *db);

//  Add a distribution to itself x-1 times ("scalar multiplication" by x)
dist_t *dist_scmul(dist_t *dr, const dist_t *da, uint64_t x);

//  multiply a distribution by a constant, spreading it by factor c
dist_t *dist_spread(dist_t *dr, const dist_t *da, int64_t c);

//  return tail -- probability mass outside (lo, hi)
double dist_tail(dist_t *d, int64_t lo, int64_t hi);

//  _FP_DIST_H_
#endif
