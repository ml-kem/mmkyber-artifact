//  mm_sample.c
//  === Simple Samplers (XXX: the Gaussian ones are placeholders.)

#include <math.h>
#include "mm_param.h"
#include "mm_sample.h"

//  Uniform sampler [0, q-1]

#define MM_Q_SZ     ((MM_LOGQ + 7) / 8)
#define MM_Q_MASK   ((1 << MM_LOGQ) - 1)

void poly_unif(int32_t *r, sha3_t *kec)
{
    uint8_t h[4] = { 0 };
    int i;
    int32_t x;

    i = 0;
    while (i < MM_D) {
        do {
            sha3_squeeze(kec, h, MM_Q_SZ);
            x = get32u_le(h) & MM_Q_MASK;
        } while (x >= MM_Q);
        r[i++] = x;
    }
}

//  Uniform sampler [-nu, nu]
/*
#define MM_NU_S     (2 * (MM_NU) + 1)
#define MM_NU_MASK  ((1 << MM_NU_BITS) - 1)

void poly_nu(int32_t *r, sha3_t *kec)
{
    uint8_t h[8] = { 0 };
    int i, l;
    int32_t x;
    uint64_t t;

    i = 0;
    l = 0;
    while (i < MM_D) {
        if (l < MM_NU_BITS) {
            sha3_squeeze(kec, h, 8);
            t = get64u_le(h);
            l = 64;
        }
        x = t & MM_NU_MASK;
        t >>= MM_NU_BITS;
        l -= MM_NU_BITS;
        if (x < MM_NU_S) {
            r[i++] = x - MM_NU;
        }
    }
}
*/

#if (MM_NU_BAR == 2)

//  Sample a polynomial in binary uniform set {0,1}.

void sample_nu(uint8_t *r, sha3_t *kec)
{
    sha3_squeeze(kec, r, MM_D / 8);
}

//  Decode a binary polynomial from bytes..

void poly_nu(int32_t *r, const uint8_t *s)
{
    int i, j, x;

    for (i = 0; i < MM_D; i += 8) {
        x = *s++;
        for (j = 0; j < 8; j++) {
            r[j] = (x >> j) & 1;
        }
        r += 8;
    }
}

#elif (MM_NU_BAR == 3)

//  Sample a polynomial in ternary uniform set {-1,0,1}.

void sample_nu(uint8_t *r, sha3_t *kec)
{
    int i;

    i = 0;
    while (i < (MM_D + 4) / 5) {
        sha3_squeeze(kec, &r[i], 1);
        if (r[i] < 243) {
            i++;
        }
    }
}

//  Decode a ternary polynomial from 52 bytes in [0, 243).

void poly_nu(int32_t *r, const uint8_t *s)
{
    int i;
    uint32_t x, a, b;

    for (i = 0; i < 51; i++) {
        x = s[i];

        a = x * 0xAAAB;
        b = (a >> 13) & 3;
        r[0] = (int32_t) b - 1;
        x = (x - b) * 0xAAAB;

        a = x * 0xAAAB;
        b = (a >> 13) & 3;
        r[1] = (int32_t) b - 1;
        x = (x - b) * 0xAAAB;

        a = x * 0xAAAB;
        b = (a >> 13) & 3;
        r[2] = (int32_t) b - 1;
        x = (x - b) * 0xAAAB;

        a = x * 0xAAAB;
        b = (a >> 13) & 3;
        r[3] = (int32_t) b - 1;
        x = (x - b) * 0xAAAB;

        a = x * 0xAAAB;
        b = (a >> 13) & 3;
        r[4] = (int32_t) b - 1;

        r += 5;
    }
    x = s[51];
    a = x * 0xAAAB;
    b = (a >> 13) & 3;
    r[0] = (int32_t) b - 1;
}

#endif

/*
=== IMPORANT NOTICE: This sampler is approximate and not constant-time;
    it is a placeholder implementation. A more appropriate Discrete Gaussian
    sampler is required in practice.
*/

//  Gaussian sampler

void poly_gauss(int32_t *r, sha3_t *kec, double sigma)
{
    int i;
    double cs2, d63;
    double x, y, w;
    uint8_t h[16] = { 0 };

    cs2 = 1.0 / 6.0 - 2.0 * (sigma * sigma);
    d63 = ldexp(1.0, -63);

    i = 0;
    while (i < MM_D) {
        sha3_squeeze(kec, h, 16);
        x = d63 * ((double) get64u_le(h)) - 1.0;
        y = d63 * ((double) get64u_le(h + 8)) - 1.0;
        w = x*x + y*y;
        if (w > 0.0 && w <= 1.0) {
            w = sqrt( cs2 * log(w) / w );
            r[i++] = rint( x * w );
            r[i++] = rint( y * w );
        }
    }
}

