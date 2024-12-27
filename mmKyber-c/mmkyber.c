//  mmkyber.c
//  === mmKyber-KEM and mmKyber-PKE implemetation

#include <string.h>

#include "mmkyber.h"
#include "mm_ring.h"
#include "mm_serial.h"
#include "mm_sample.h"

#include "sha3_t.h"

//  use NTT-less decryption/decapsulation
//#define MM_NO_NTT_DEC

//  for indexing the A matrix
#define MM_A_IDX(i,j) (((i) * MM_N + (j)) * MM_D)

//  helper function for setting up SHAKE for output

static void kec_setup(sha3_t *kec, const uint8_t *buf, size_t len)
{
    //  use SHAKE128 for security level 128 and below, otherwise SHAKE256
#ifdef MM_128
    sha3_init(kec, SHAKE128_RATE);
#else
    sha3_init(kec, SHAKE256_RATE);
#endif

    sha3_absorb(kec, buf, len);
    sha3_pad(kec, SHAKE_PAD);
}

//  mmKEM & mmPKE: mmSetup(1^lambda, N): Generate public parameter A from seed.

void mm_setup(int32_t *a, const uint8_t seed_a[16])
{
    int i, j;
    uint8_t seed[24];
    sha3_t kec;

    memcpy(seed, seed_a, 16);

    for (i = 0; i < MM_M; i++) {
        for (j = 0; j < MM_N; j++) {

            seed[16] = i;
            seed[17] = j;
            seed[18] = 'A';

            //  A matrix generation is always SHAKE128
            sha3_init(&kec, SHAKE128_RATE);
            sha3_absorb(&kec, seed, 19);
            sha3_pad(&kec, SHAKE_PAD);
            poly_unif(&a[MM_A_IDX(i,j)], &kec);
        }
    }
}

//  mmKEM & mmPKE: mmKGen(pp): Generate individual public key.

size_t mm_kgen( uint8_t *pk, uint8_t *sk,
                const int32_t *a_mat, const uint8_t seed_k[32])
{
    int i, j;
    int32_t s[MM_M][MM_D];
    int32_t e[MM_N][MM_D];
    int32_t b[MM_D];
    size_t pk_sz;
    uint8_t seed[40];
    uint8_t buf[MM_NU_SZ];
    sha3_t kec;

    //  (s, e) <- U(Snu^m) x U(Snu^n)
    memcpy(seed, seed_k, 32);
    seed[32] = 'S';
    kec_setup(&kec, seed, 33);
    for (i = 0; i < MM_M; i++) {
        sample_nu(sk, &kec);
        poly_nu(s[i], sk);
        sk += MM_NU_SZ;
        polyr_fntt(s[i]);
    }

    //  b := A^T * s + e
    seed[32] = 'E';
    kec_setup(&kec, seed, 33);
    for (i = 0; i < MM_N; i++) {
        sample_nu(buf, &kec);
        poly_nu(e[i], buf);
        polyr_fntt(e[i]);
    }

    pk_sz = 0;
    for (i = 0; i < MM_N; i++) {

        polyr_zero(b);
        for (j = 0; j < MM_M; j++) {
            polyr_ntt_mul_add(b, &a_mat[MM_A_IDX(j, i)], s[j]);
        }

        polyr_scale(b, MONT_RR, b);     //  remove the montgomery factor
        polyr_add(b, b, e[i]);
        polyr_norm(b);

        //  t := b
        pk_sz += poly_serial(pk + pk_sz, b, MM_LOGQ);
    }

    return pk_sz;
}


//  mmKEM & mmPKE private: mmEnc^i(pp; r): Shared ciphertext.

static size_t mm_enc_i( uint8_t *ct, const int32_t *a_mat,
                        const int32_t r[][MM_D],        //  ntt domain
                        const int32_t e[][MM_D])        //  normal domain
{
    int i, j;
    int32_t c[MM_D];                //  no need to store the whole vector
    size_t ct_sz;

    //  c := A * r + e_u
    ct_sz = 0;
    for (i = 0; i < MM_M; i++) {
        polyr_zero(c);
        for (j = 0; j < MM_N; j++) {
            polyr_ntt_mul_add(c, &a_mat[MM_A_IDX(i, j)], r[j]);
        }
        polyr_intt(c);
        polyr_add(c, c, e[i]);

        //  u := [c mod q] (2^du)
        ct_sz += poly_compress(ct + ct_sz, c, MM_DU);
    }

    return ct_sz;
}

//  mmKEM private: mmEncap^d(pp, pk_i; r, r_i): Individual ciphertext.

static size_t mm_encap_d(   uint8_t *ct, uint8_t *k,
                            const uint8_t *pk,
                            const int32_t r[][MM_D], // in ntt domain
                            const int32_t y[MM_D ] )
{
    int i;
    int32_t b[MM_D];
    int32_t c[MM_D];

    polyr_zero(c);
    for (i = 0; i < MM_N; i++) {

        //  b'_i := t_i
        pk += poly_deserial(b, pk, MM_LOGQ);

        //  c_i := < b'_i, r > + y_i
        polyr_ntt_mul_add(c, b, r[i]);
    }
    polyr_intt(c);
    polyr_add(c, c, y);

    //  fast generation
    poly_gen_ct_k(ct, k, c);

    return MMKEM_CTI_SZ;
}

//  mmKEM: mmEncap(pp, (pk_i) for i in [N]): Encapsulate to N recipients.

size_t mm_encap(uint8_t *ct, uint8_t *kk,
                const int32_t *a_mat, const uint8_t *pk[],
                const uint8_t seed_e[32], size_t n)
{
    size_t i;
    sha3_t kec;
    int32_t r_u[MM_N][MM_D];
    int32_t e_u[MM_M][MM_D];

    uint8_t buf[64];
    size_t ct_sz;

    //  r := (r, e_u) <= D^n_sigma0 x D^M_sigma0
    memcpy(buf, seed_e, 32);

    buf[32] = 'R';
    kec_setup(&kec, buf, 33);
    for (i = 0; i < MM_N; i++) {
        poly_gauss(r_u[i], &kec, MM_SIGMA0);
        polyr_fntt(r_u[i]);
    }

    buf[32] = 'e';
    kec_setup(&kec, buf, 33);
    for (i = 0; i < MM_M; i++) {
        poly_gauss(e_u[i], &kec, MM_SIGMA0);
    }

    //  ^ct <- mmEnc^i(pp; r)
    ct_sz = mm_enc_i(ct, a_mat, r_u, e_u);

    int32_t y[MM_D];

    for (i = 0; i < n; i++) {

        put64u_le(buf + 32, i);
        buf[40] = 'r';
        kec_setup(&kec, buf, 41);
        poly_gauss(y, &kec, MM_SIGMA1);

        //  (~ct_i, K_i) <- mmEncap^d(pp. pk_i; r, r_i)
        mm_encap_d(ct + ct_sz, kk, /* a_mat, */ pk[i], r_u, y /*, &kec */);
        kk += MMKEM_K_SZ;
        ct_sz += MMKEM_CTI_SZ;
    }

    return ct_sz;
}

//  mmKEM: mmDecap(pp, sk, ct): Decapsulate individual ciphertext (ctu,cti).

#ifdef MM_NO_NTT_DEC

//  No-NTT multiply with a ternary secret.

static void poly_mul1_add16(uint16_t *r, uint16_t *f, int32_t *g)
{
    int i, j, x;

    for (i = 0; i < MM_D; i++) {
        x = g[i];
        for (j = 0; j < MM_D; j++) {
            r[j] += x * f[j];
        }
        r++;
    }
}

void mm_decap(  uint8_t *k, const uint8_t *sk,
                const uint8_t *ctu, const uint8_t *cti)
{
    int i;
    int32_t s[MM_D];
    uint16_t c[MM_D], w[2 * MM_D];
    uint16_t x, b;

    memset(w, 0, sizeof(w));
    for (i = 0; i < MM_M; i++) {

        //  Expand private key
        poly_nu(s, sk);
        sk += MM_NU_SZ;

        //  c' := u mod 2^du
        ctu += poly_deserial16(c, ctu, MM_DU, MM_D);

        //  w := <c',s> mod 2^u_i
        poly_mul1_add16(w, c, s);
    }

    memset(k, 0, MMKEM_K_SZ);
    for (i = 0; i < MMKEM_K_SZ * 8; i++) {
        b = (cti[i >> 3] >> (i & 7)) & 1;
        x = (w[i] - w[i + MM_D]) >> (MM_DU - 3);
        x = ((x + 2 * b + 1) >> 2) & 1;
        k[i >> 3] |= x << (i & 7);
    }
}

#else

//  Uses NTT (the partial sums fit into q)

void mm_decap(  uint8_t *k, const uint8_t *sk,
                const uint8_t *ctu, const uint8_t *cti)
{
    int i;
    int32_t x, s[MM_D], c[MM_D], w[MM_D];
    uint16_t b;

    polyr_zero(w);
    for (i = 0; i < MM_M; i++) {

        //  Expand private key
        poly_nu(s, sk);
        sk += MM_NU_SZ;
        polyr_fntt(s);

        //  c' := u mod 2^du
        ctu += poly_deserial(c, ctu, MM_DU);
        polyr_fntt(c);

        //  w := <c',s> mod 2^u_i
        polyr_ntt_mul_add(w, c, s);
    }
    polyr_intt(w);

    memset(k, 0, MMKEM_K_SZ);
    for (i = 0; i < MMKEM_K_SZ * 8; i++) {

        //  make <c',s> signed and clip to 3-bit range
        x   = w[i];
        x   -=  ~((x - (MM_Q / 2)) >> 31) & MM_Q;
        x   >>= MM_DU - 3;

        //  rec
        b = (cti[i >> 3] >> (i & 7)) & 1;
        x = ((x + 2 * b + 1) >> 2) & 1;
        k[i >> 3] |= x << (i & 7);
    }
}

#endif


//  mmPKE private: mmEnc^d(pp, pk_i, m_i; r, r_i)

static size_t mm_enc_d( uint8_t *ct,
                        const uint8_t *pk,
                        const uint8_t *m,
                        const int32_t r[][MM_D], // in ntt domain
                        const int32_t y[MM_D])
{
    int i, j, x;
    int32_t b[MM_D];
    int32_t c[MM_D];

    polyr_zero(c);
    for (i = 0; i < MM_N; i++) {

        //  b'_i := t_i
        pk += poly_deserial(b, pk, MM_LOGQ);

        //  c_i := < b'_i, r > + y_i
        polyr_ntt_mul_add(c, b, r[i]);
    }
    polyr_intt(c);
    polyr_add(c, c, y);

    //  encode message
    for (i = 0; i < MMPKE_M_SZ; i++) {
        x = m[i];
        for (j = 0; j < 8; j++) {
            //  [q/2]*m_i
            c[8 * i + j] += (-((x >> j) & 1)) & ((MM_Q + 1) / 2);
        }
    }
    //  normalize
    polyr_norm(c);

    //  compress to dv bits, return length
    return poly_compress(ct, c, MMPKE_DV);
}


//  mmPKE: mmEnc(pp, (pk_i), (m_i) for i in [N]): Encrypt to N recipients.

size_t mm_enc(  uint8_t *ct, const int32_t *a_mat,
                const uint8_t *pk[], const uint8_t *mm,
                const uint8_t seed_e[32], size_t n)
{
    size_t i;
    sha3_t kec;
    int32_t r_u[MM_N][MM_D];
    int32_t e_u[MM_M][MM_D];

    uint8_t buf[64];
    size_t ct_sz;

    //  r := (r, e_u) <= D^n_sigma0 x D^M_sigma0
    memcpy(buf, seed_e, 32);

    buf[32] = 'R';
    kec_setup(&kec, buf, 33);
    for (i = 0; i < MM_N; i++) {
        poly_gauss(r_u[i], &kec, MM_SIGMA0);
        polyr_fntt(r_u[i]);
    }

    buf[32] = 'e';
    kec_setup(&kec, buf, 33);
    for (i = 0; i < MM_M; i++) {
        poly_gauss(e_u[i], &kec, MM_SIGMA0);
    }

    //  ^ct <- mmEnc^i(pp; r)
    ct_sz = mm_enc_i(ct, a_mat, r_u, e_u);

    int32_t y[MM_D];

    for (i = 0; i < n; i++) {

        //  r_i := y_i <- D_sigma1
        put64u_le(buf + 32, i);
        buf[40] = 'r';
        kec_setup(&kec, buf, 41);
        poly_gauss(y, &kec, MM_SIGMA1);

        //  (~ct_i, K_i) <- mmEncap^d(pp. pk_i; r, r_i)
        mm_enc_d(ct + ct_sz, pk[i], mm, r_u, y);
        ct_sz += MMPKE_CTI_SZ;
        mm  += MMPKE_M_SZ;
    }

    return ct_sz;
}

//  mmPKE: mmDec(pp, sk, ct): Decrypt a message

#ifdef MM_NO_NTT_DEC

//  NTT-free version

void mm_dec(uint8_t *m, const uint8_t *sk,
            const uint8_t *ctu, const uint8_t *cti)
{
    int i, x;
    int32_t s[MM_D];
    uint16_t u[MM_D], v[MM_D], w[2 * MM_D];

    //  u' := [u mod 2^dv](2^du)
    poly_deserial16(v, cti, MMPKE_DV, MMPKE_M_SZ * 8);

    memset(w, 0, sizeof(w));
    for (i = 0; i < MM_M; i++) {

        //  Expand private key
        poly_nu(s, sk);
        sk += MM_NU_SZ;

        //  c' := [u mod 2^du](q)
        ctu += poly_deserial16(u, ctu, MM_DU, MM_D);

        //  w := <u,s> mod 2^u_i
        poly_mul1_add16(w, u, s);
    }

    //  m := [ u' - <u,s> mod 2^d_u ]_2
    memset(m, 0, MMPKE_M_SZ);
    for (i = 0; i < MMPKE_M_SZ * 8; i++) {
        x = (v[i] << (MM_DU - MMPKE_DV)) - (w[i] - w[i + MM_D]);
        x = ((x + (1 << (MM_DU - 2))) >> (MM_DU - 1)) & 1;
        m[i >> 3] |= (x & 1) << (i & 7);
    }
}

#else

//  This version uses NTT

void mm_dec(uint8_t *m, const uint8_t *sk,
            const uint8_t *ctu, const uint8_t *cti)
{
    int i, x;
    int32_t s[MM_D], u[MM_D], w[MM_D];
    uint16_t v[MM_D];

    //  u' := [u mod 2^dv](2^du)
    poly_deserial16(v, cti, MMPKE_DV, MMPKE_M_SZ * 8);

    polyr_zero(w);
    for (i = 0; i < MM_M; i++) {

        //  Expand private key
        poly_nu(s, sk);
        sk += MM_NU_SZ;
        polyr_fntt(s);

        //  u' := [u mod 2^du](q)
        ctu += poly_deserial(u, ctu, MM_DU);
        polyr_fntt(u);

        //  w := <u,s> mod 2^u_i
        polyr_ntt_mul_add(w, u, s);
    }
    polyr_intt(w);

    memset(m, 0, MMPKE_M_SZ);
    for (i = 0; i < MMPKE_M_SZ * 8; i++) {

        x   = w[i];
        //  make <u,s> signed and clip to 2**du range
        x   -=  ~((x - (MM_Q / 2)) >> 31) & MM_Q;
        x   &= (1 << MM_DU) - 1;

        //  m := [ u' - <u,s> mod 2^d_u ]_2
        x = (v[i] << (MM_DU - MMPKE_DV)) - x;
        x = ((x + (1 << (MM_DU - 2))) >> (MM_DU - 1)) & 1;
        m[i >> 3] |= (x & 1) << (i & 7);
    }
}

#endif
