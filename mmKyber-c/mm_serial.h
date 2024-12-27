//  mm_serial.h
//  === Bit packing, compression, and Ding-Peikert Reconciliation.

#ifndef _MM_SERIAL_H_
#define _MM_SERIAL_H_

#include "plat_local.h"
#include "mm_param.h"

//  Pack dx-bit elements from "p" to "b". Return byte length written to "b".

static inline
size_t poly_serial(uint8_t *b, const int32_t *p, int dx)
{
    int64_t x, t, m;
    int i, j, l;

    m = (1 << dx) - 1;
    l = 0;
    t = 0;
    j = 0;
    for (i = 0; i < MM_D; i++) {

        //  read coefficint
        x = p[i] & m;

        //  serialize
        t |= x << l;
        l += dx;
        while (l >= 8) {
            b[j++] = t & 0xFF;
            t >>= 8;
            l -= 8;
        }
    }
    //  remaining part?
    if (l > 0) {
        b[j++] = t & 0xFF;
    }

    return (size_t) j;
}


//  Unpack dx-bit elements from "b" to "p". Return bytes read from "b".

static inline
size_t poly_deserial(int32_t *p, const uint8_t *b, int dx)
{
    uint32_t t, m, s;
    int i, j, l;

    if (dx == MM_LOGQ) {        //  signedness
        s = 0;
    } else {
        s = 1 << (dx - 1);
    }

    m = (1 << dx) - 1;
    i = 0;
    l = 0;
    t = 0;
    for (j = 0; j < MM_D; j++) {

        //  decode a coefficient
        while (l < dx) {
            t |= ((uint32_t) b[i++]) << l;
            l += 8;
        }
        p[j] = ((t + s) & m) - s;
        l -= dx;
        t >>= dx;
    }

    return (size_t) i;
}


//  Unpack "len" dx-bit elements from "b" to "p" where p is 16-bits wide.

static inline
size_t poly_deserial16(uint16_t *p, const uint8_t *b, int dx, int len)
{
    uint32_t t, m;
    int i, j, l;

    m = (1 << dx) - 1;
    i = 0;
    l = 0;
    t = 0;
    for (j = 0; j < len; j++) {

        //  decode a coefficient
        while (l < dx) {
            t |= ((uint32_t) b[i++]) << l;
            l += 8;
        }
        p[j] = t & m;
        l -= dx;
        t >>= dx;
    }

    return (size_t) i;
}


//  Compress a polynomial "a" to bytes in "b", "dx" bits per coefficient.
//  Return number of bytes written to "b".

static inline
size_t poly_compress(uint8_t *b, const int32_t *p, int dx)
{
    int64_t d, x, y;
    int64_t t, m;
    int i, j, l;

    m = (1 << dx) - 1;
    l = 0;
    t = 0;
    j = 0;
    for (i = 0; i < MM_D; i++) {

        //  scale, add a rounding constant
        x = (((int64_t) p[i]) << dx) + (MM_Q / 2l);

        //  divide by Q, constant time (x < 2**51)
        d = x >> MM_LOGQ;
        y = x - d * MM_Q;
        d += y >> MM_LOGQ;
        y = x - d * MM_Q;
        d += ((y - MM_Q) >> 31) + 1;

        //  serialize
        t |= (d & m) << l;
        l += dx;
        while (l >= 8) {
            b[j++] = t & 0xFF;
            t >>= 8;
            l -= 8;
        }
    }
    //  remaining part?
    if (l > 0) {
        b[j++] = t & 0xFF;
    }

    return (size_t) j;
}


//  Decompress a polynomial from "b" ("dx" bits per coefficient) into "p",
//  Return number of bytes read from "b".

static inline
size_t poly_decompress(int32_t *p, const uint8_t *b, int dx)
{
    int64_t x, t, m;
    int i, j, l;

    m = (1 << dx) - 1;
    i = 0;
    l = 0;
    t = 0;
    for (j = 0; j < MM_D; j++) {

        //  decode a coefficient
        while (l < dx) {
            t |= ((int64_t) b[i++]) << l;
            l += 8;
        }
        x = t & m;
        l -= dx;
        t >>= dx;

        //  round up and store coefficient
        x = (x * MM_Q) + (1l << (dx - 1));
        p[j] = (int32_t) (x >> dx);
    }

    return (size_t) i;
}

//  Create ciphertext and key bits from "approximate shared secret."

static inline
void poly_gen_ct_k(uint8_t *ct, uint8_t *k, int32_t *c)
{
    int i, j;
    int32_t x, d;
    uint8_t a, b;

    for (i = 0; i < MMKEM_CTI_SZ; i++) {
        a = 0;
        b = 0;
        for (j = 0; j < 8; j++) {

            x = 4 * c[8 * i + j];

            //  divide by q -- constant time
            d = x >> MM_LOGQ;
            x -= d * MM_Q;
            d += ((x - MM_Q) >> 31) + 1;

            a |= (d & 1) << j;
            b |= (((d + 1) >> 1) & 1) << j;
        }
        ct[i] = a;
        k[i] = b;
    }
}

#endif
