//  mm_ring.h
//  === Header: General Polynomial Ring Arithmetic.

#ifndef _MM_RING_H_
#define _MM_RING_H_

#include "plat_local.h"
#include "mm_param.h"

/*
d = 256
q = 2^25 - 2^12 + 1
r = 2^32 % q
rr = r^2 % q
di = lift(rr * Mod(d,q)^-1)
qi = lift(Mod(-q,2^32)^-1)
*/

//  derived parameters
#define MONT_Q      MM_Q
#define MONT_LOGQ   MM_LOGQ
#define MONT_R      524160l
#define MONT_RR     33546244l
#define MONT_DI     33157153l
#define MONT_QI     16773119l

//  fast normalizing reduction
static inline int32_t mont_red1(int32_t x)
{
    //  (no input restrictions)
    x -= (x >> MONT_LOGQ) * MONT_Q;

    //  these are the exact bounds for 32-bit signed x
    XASSERT(x >= -0x3FFC0);     //  -262079
    XASSERT(x <= 0x203EFC0);    //  q + 262079
    return x;
}

//  Montgomery reduction. Returns r in [-q,q-1] so that r == (x/2^32) mod q.

static inline int32_t mont_redc(int64_t x)
{
    int32_t r;

    //  prove these input bounds; -(8*q)^2 < x < (8*q)^2
    XASSUME(x >= -(64l * MONT_Q * MONT_Q));
    XASSUME(x <= (64l * MONT_Q * MONT_Q));

    r = (int32_t) x * MONT_QI;
    r = (x + ((int64_t) r) * ((int64_t) MONT_Q)) >> 32;

    //  prove output bounds (only one coditional addition is required)
    XASSERT(r >= -MONT_Q);
    XASSERT(r < MONT_Q);

    return r;
}

//  Montgomery multiplication

static inline int32_t mont_mulq(int32_t x, int32_t y)
{
    int32_t r;

    XASSUME(x >= -8 * MONT_Q);
    XASSUME(x <= 8 * MONT_Q);
    XASSUME(y >= -8 * MONT_Q);
    XASSUME(y <= 8 * MONT_Q);

    r = mont_redc(((int64_t) x) * ((int64_t) y));

    XASSERT(r >= -MONT_Q);
    XASSERT(r < MONT_Q);

    return r;
}

//  Add q conditionally if negative

static inline int32_t mont_cadd(int64_t x)
{
    int32_t r;

    XASSUME(x >= -MONT_Q);
    XASSUME(x < MONT_Q);

    r = x + (MONT_Q & (x >> 31));

    XASSERT(r >= 0);
    XASSERT(r < MONT_Q);
    XASSERT(r == x || r == x + MONT_Q);

    return r;
}

//  === polyr.c prototypes

//  Zero polynomial
static inline void polyr_zero(int32_t *r)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        r[i] = 0;
    }
}

//  Copy polynomial
static inline void polyr_copy(int32_t *r, const int32_t *a)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        r[i] = a[i];
    }
}

//  Add polynomials:  r = f + g
static inline void polyr_add(int32_t *r, const int32_t *f, const int32_t *g)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        r[i] = f[i] + g[i];
    }
}

//  Subtract polynomials:  r = f - g
static inline void polyr_sub(int32_t *r, const int32_t *f, const int32_t *g)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        r[i] = f[i] - g[i];
    }
}

//  Single conditional add, puts [-q,q] -> [0,q]
static inline void polyr_cadd(int32_t *r)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        r[i] = mont_cadd(r[i]);
    }
}

//  transform into range [0,q-1]

static inline void polyr_norm(int32_t *r)
{
    int i;
    int32_t x;

    for (i = 0; i < MM_D; i++) {
        x = mont_red1(r[i]);
        x = mont_cadd(x);
        x = mont_cadd(x - MONT_Q);
        r[i] = x;
    }
}

//  scalar multiplication c * f

static inline void polyr_scale(int32_t *r, int32_t c, const int32_t *f)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        r[i] = mont_mulq(f[i], c);
    }
}

//  Coefficient multiply:  fg = f * g,  Montgomery reduction.

static inline void polyr_ntt_mul(   int32_t *fg,
                                    const int32_t *f, const int32_t *g)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        fg[i] = mont_mulq(f[i], g[i]);
    }
}

//  Coefficient multiply and add:  r += a * b, Montgomery reduction.

static inline void  polyr_ntt_mul_add(  int32_t *fg,
                                        const int32_t *f, const int32_t *g)
{
    int i;

    for (i = 0; i < MM_D; i++) {
        fg[i] += mont_mulq(f[i], g[i]);
    }
}

//  === Protytpes -- polyr.c

//  Forward NTT (negacyclic -- evaluate polynomial at factors of x^n+1).
void polyr_fntt(int32_t *f);

//  Reverse NTT (negacyclic -- x^n+1), normalize by 1/(n*r).
void polyr_intt(int32_t *f);

#endif
