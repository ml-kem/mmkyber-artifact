//  crt_ntt.c
//  === NTT Convolutions

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

#include "crt_ntt.h"

#define CRT_NTT_PRIME_BASE 32

//  my types
typedef uint64_t plimb_t;
typedef __uint128_t plimb2_t;

//  prime base -- q = 2^64 - 2^32 * i + 1

static const plimb_t crt_ntt_q[CRT_NTT_PRIME_BASE] = {
    0xffffffff00000001llu, // [ 0] i=   1
    0xfffffffc00000001llu, // [ 1] i=   4
    0xffffffd300000001llu, // [ 2] i=  45
    0xffffffca00000001llu, // [ 3] i=  54
    0xffffffc600000001llu, // [ 4] i=  58
    0xffffffb500000001llu, // [ 5] i=  75
    0xffffffb200000001llu, // [ 6] i=  78
    0xffffffa300000001llu, // [ 7] i=  93
    0xffffff9300000001llu, // [ 8] i= 109
    0xffffff8400000001llu, // [ 9] i= 124
    0xffffff8200000001llu, // [10] i= 126
    0xffffff7300000001llu, // [11] i= 141
    0xffffff6400000001llu, // [12] i= 156
    0xffffff5700000001llu, // [13] i= 169
    0xffffff3700000001llu, // [14] i= 201
    0xffffff0f00000001llu, // [15] i= 241
    0xffffff0000000001llu, // [16] i= 256
    0xfffffefb00000001llu, // [17] i= 261
    0xfffffeef00000001llu, // [18] i= 273
    0xfffffee500000001llu, // [19] i= 283
    0xfffffeda00000001llu, // [20] i= 294
    0xfffffed100000001llu, // [21] i= 303
    0xfffffe9a00000001llu, // [22] i= 358
    0xfffffe2e00000001llu, // [23] i= 466
    0xfffffe1d00000001llu, // [24] i= 483
    0xfffffdea00000001llu, // [25] i= 534
    0xfffffde000000001llu, // [26] i= 544
    0xfffffdd400000001llu, // [27] i= 556
    0xfffffdb400000001llu, // [28] i= 588
    0xfffffdb100000001llu, // [29] i= 591
    0xfffffdad00000001llu, // [30] i= 595
    0xfffffda500000001llu, // [31] i= 603
};

//  forward and inverse transform twiddle values

static plimb_t crt_ntt_cf[CRT_NTT_PRIME_BASE][32];
static plimb_t crt_ntt_ci[CRT_NTT_PRIME_BASE][32];

//  Basic arithmetic primitives
//  Written as separate functions -- so that they can be instrumented easily.

static inline plimb_t crt_ntt_addq(plimb_t x, plimb_t y, plimb_t q)
{
    plimb2_t z = ((plimb2_t) x) + ((plimb2_t) y);
    return z % ((plimb2_t) q);
}

static inline plimb_t crt_ntt_subq(plimb_t x, plimb_t y, plimb_t q)
{
    plimb2_t z = ((plimb2_t) q) + ((plimb2_t) x) - ((plimb2_t) y);
    return z % ((plimb2_t) q);
}

static inline plimb_t crt_ntt_mulq(plimb_t x, plimb_t y, plimb_t q)
{
    return (((plimb2_t) x) * ((plimb2_t) y)) % ((plimb2_t) q);
}

static inline plimb_t crt_ntt_sqrq(plimb_t x, plimb_t q)
{
    return (((plimb2_t) x) * ((plimb2_t) x)) % ((plimb2_t)q);
}

//  Compute x**e (mod q).

static plimb_t crt_ntt_expq(plimb_t x, plimb_t e, plimb_t q)
{
    plimb_t y;

    y = 1;
    if (e & 1) {
        y = x;
    }
    e >>= 1;

    while (e > 0) {
        x = crt_ntt_sqrq(x, q);
        if (e & 1)
            y = crt_ntt_mulq(x, y, q);
        e >>= 1;
    }

    return y;
}

//  test multiplicative order == 0

static plimb_t crt_ntt_ord(plimb_t g, size_t n, plimb_t q)
{
    size_t i;
    plimb_t x, y;

    //  reduce order to 2n
    x = crt_ntt_expq(g, (q - 1) / (n << 1), q);

    y = x;  //  test order via repeated squaring
    for (i = 1; i < 2 * n; i <<= 1) {
        y = crt_ntt_sqrq(y, q);
        if (y == 1)
            break;
    }
    if (i != n)
        return 0;  //   return 0 on failure

    return x;
}


//  Forward NTT (cyclic)

static void crt_ntt_fntt(   plimb_t *f, const plimb_t *cf,
                            size_t n, plimb_t q)
{
    size_t i, j, k;
    plimb_t x, y, z;
    plimb_t *p0, *p1, *p2;

    for (k = 1, j = n >> 1; j > 0; k <<= 1, j >>= 1) {

        p0 = f;
        p1 = p0 + j;
        p2 = p1 + j;

        while (p1 < p2) {
            x = *p0;
            y = *p1;
            *p0++ = crt_ntt_addq(x, y, q);
            *p1++ = crt_ntt_subq(x, y, q);
        }
        p0 = p2;
        p1 = p0 + j;
        p2 = p1 + j;

        z = 1;
        for (i = 1; i < k; i++) {

            if (i & 1) {
                z = crt_ntt_mulq(z, cf[0], q);
            } else {
                z = crt_ntt_mulq(z, cf[ __builtin_ctzll(i) ], q);
            }

            while (p1 < p2) {
                x = *p0;
                y = *p1;
                y = crt_ntt_mulq(y, z, q);
                *p0++ = crt_ntt_addq(x, y, q);
                *p1++ = crt_ntt_subq(x, y, q);
            }
            p0 = p2;
            p1 = p0 + j;
            p2 = p1 + j;
        }
    }
}

//  Reverse NTT (cyclic)

static void crt_ntt_intt(   plimb_t *f, const plimb_t *cr,
                            size_t n, plimb_t q)
{
    size_t i, j, k;
    plimb_t x, y, z;
    plimb_t *p0, *p1, *p2;

    for (j = 1, k = n >> 1; k > 0; j <<= 1, k >>= 1) {

        p0 = f;
        p1 = p0 + j;
        p2 = p1 + j;

        while (p1 < p2) {
            x = *p0;
            y = *p1;
            *p0++ = crt_ntt_addq(x, y, q);
            *p1++ = crt_ntt_subq(x, y, q);
        }

        z = q - 1;
        for (i = 1; i < k; i++) {

            if (i & 1) {
                z = crt_ntt_mulq(z, cr[0], q);
            } else {
                z = crt_ntt_mulq(z, cr[ __builtin_ctzll(i) ], q);
            }

            p0 = p2;
            p1 = p0 + j;
            p2 = p1 + j;

            while (p1 < p2) {
                x = *p0;
                y = *p1;
                *p0++ = crt_ntt_addq(x, y, q);
                y     = crt_ntt_subq(y, x, q);
                *p1++ = crt_ntt_mulq(y, z, q);
            }
        }
    }
}

//  Elementvise vector product  r = f (*) g.

static void crt_ntt_xmul(   plimb_t *fg,
                            const plimb_t *f, const plimb_t *g,
                            size_t n, plimb_t q)
{
    size_t i;

    //  multiply each element point-by-point
    for (i = 0; i < n; i++) {
        fg[i] = crt_ntt_mulq(f[i], g[i], q);
    }
}

/*
//  Multiply with a scalar  cf = c * f.

static void crt_ntt_cmul(   plimb_t *cf, plimb_t c, const plimb_t *f,
                            size_t n, plimb_t q)
{
    size_t i;

    for (i = 0; i < n; i++) {
        cf[i] = crt_ntt_mulq(c, f[i], q);
    }
}
*/

//  build the log-sized twiddle table

static void pow3_tab(   plimb_t *lv, plimb_t h, size_t l, plimb_t q)
{
    plimb_t h3;

    h3 = crt_ntt_mulq( crt_ntt_sqrq(h, q), h, q );  //  h3 = h^3
    while (l > 1) {
        lv[--l] = q - h3;
        h3 = crt_ntt_sqrq( h3, q );
    }
    lv[0] = q - h3;
}

//  initialize tables

size_t crt_ntt_init()
{
    size_t i;
    plimb_t q = 1, h = 0, hi = 0, g = 0;
    const size_t l = 31;

    for (i = 0; i < CRT_NTT_PRIME_BASE; i++) {

        //  find generator
        q = crt_ntt_q[i];
        h = 0;
        for (g = 2; g < q; g++) {
            h = crt_ntt_ord(g, 1lu << l, q);
            if (h > 1)
                break;
        }
        if (g == q) {
            printf("No %lu-order generator found for %lu\n", 2lu << l, q);
            exit(0);
        }
        hi = crt_ntt_expq(h, q - 2, q); //  hi = 1/h

        //  printf("q= %20lu; h= %20lu;\n", q, h);

        //  make tables
        pow3_tab(crt_ntt_cf[i], h, l, q);   //  log-sized twiddle table, forward
        pow3_tab(crt_ntt_ci[i], hi, l, q); //   log-sized twiddle table, inverse
    }

    return 64 * CRT_NTT_PRIME_BASE - 1;
}

//  currently does nothing

void crt_ntt_free()
{
    ;
}

//  compute a convolution r += f * g. r needs to hold "bits"

size_t crt_ntt_convol(  mpz_t r[],
                        mpz_t f[], size_t f_len,
                        mpz_t g[], size_t g_len,
                        size_t bits)
{
    size_t  i, j, n, b, r_len;
    plimb_t *fv, *gv, x, q;
    mpz_t   bz, cz, xz;

    mpz_inits(bz, cz, xz, NULL);

    assert(f_len > 0 && g_len > 0);
    r_len = f_len + g_len - 1;

    n = 1;                              //  ntt size
    while (n < r_len) {
        n <<= 1;
    }
    fv = calloc(n, sizeof(plimb_t));
    if (f == g) {
        gv = fv;
    } else {
        gv = calloc(n, sizeof(plimb_t));
    }

    //  build prime base modulus

    mpz_set_ui(bz, crt_ntt_q[0]);
    b = 1;
    while (mpz_sizeinbase(bz, 2) < bits) {
        mpz_mul_ui(bz, bz, crt_ntt_q[b++]);
    }

    //  CRT NTT's

    for (i = 0; i < b; i++) {

        q   = crt_ntt_q[i];

        //  transform f
        for (j = 0; j < f_len; j++) {
            fv[j] = mpz_fdiv_ui(f[j], q);
        }
        for (j = f_len; j < n; j++) {
            fv[j] = 0;
        }
        crt_ntt_fntt(   fv, crt_ntt_cf[i], n, q);

        //  transform g
        if (f != g) {
            for (j = 0; j < g_len; j++) {
                gv[j] = mpz_fdiv_ui(g[j], q);
            }
            for (j = g_len; j < n; j++) {
                gv[j] = 0;
            }
            crt_ntt_fntt(   gv, crt_ntt_cf[i], n, q);
        }

        //  multiply into f = f * g
        crt_ntt_xmul(   fv, fv, gv,     n, q);

        //  inverse transform
        crt_ntt_intt(   fv, crt_ntt_ci[i], n, q);

        //  build cofactor coefficient
        mpz_divexact_ui(cz, bz, q);     //  c = b / q
        x = mpz_fdiv_ui(cz, q);         //  x = c % q
        x = crt_ntt_mulq(n, x, q);          //  x = (n * c) % q
        x = crt_ntt_expq(x, q - 2, q);      //  x = 1/(n*c)  mod q
        mpz_mul_ui(cz, cz, x);          //  c = x * (b / q)

        //  add to r (mod b)
        for (j = 0; j < r_len; j++) {
            mpz_mul_ui(xz, cz, fv[j]);  //  x = c * f[j]
            mpz_add(xz, xz, r[j]);
            mpz_mod(r[j], xz, bz);      //  r[j] = (r[j] + c * x) % b
        }
    }
    if (gv != fv) {
        free(gv);
    }
    free(fv);

    mpz_clears(bz, cz, xz, NULL);

    return r_len;
}

