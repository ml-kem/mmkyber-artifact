//  test_main.c
//  === test and benchmark code

#if defined(BENCH) || defined(TESTVEC)

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include "plat_local.h"
#include "sha3_t.h"
#include "mmkyber.h"
#include "mm_param.h"

#ifndef MM_N_MAX
#define MM_N_MAX 1024
#endif

#ifndef MM_REP_TOT
#define MM_REP_TOT 10240
#endif

//  used for debug stuff

void dbg_hex(const void *b, size_t b_sz, const char *lab)
{
    size_t i;
    printf("%s[%zu] = ", lab, b_sz);
    for (i = 0; i < b_sz; i++) {
        printf("%02x", ((const uint8_t *) b)[i]);
    }
    printf("\n");
}

void dbg_vec(const int32_t *v, size_t l, const char *lab)
{
    size_t i;
    printf("%s[%zu] = {", lab, l);
    for (i = 0; i < l; i++) {
        printf(" %d,", v[i]);
    }
    printf(" };\n");
}

uint32_t dbg_sum(const void *b, size_t b_sz, const char *lab)
{
    size_t i;
    uint32_t x;

    x = 1;
    for (i = 0; i < b_sz; i++) {
        x = x * 0x103 + ((const uint8_t *) b)[i];
    }
    printf("%s[%zu] chk %08x\n", lab, b_sz, x);

    return x;
}

static double get_sec()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double) tv.tv_sec) + 1E-6*((double) tv.tv_usec);
}

int main()
{
    //  for our "test vectors"
    uint8_t seed_a[] = "0123456789abcdef";
    uint8_t seed_k[] = "000102030405060708090a0b0c0d0e0f";
    uint8_t seed_e[] = "00112233445566778899aabbccddeeff";

    uint8_t *p, *pk[MM_N_MAX], *sk[MM_N_MAX];
    uint64_t cc = 0;    //  cycle counts
    double dd = 0.0;    //  seconds
    int i, nn, rep, iter;
    size_t  nn_ct_sz = 0;
#ifdef MM_PKE
    uint8_t buf[64];
#endif

    //  we have these as global
    int32_t a_mat[MM_M * MM_N * MM_D];
    uint8_t kdata[MM_N_MAX * (MM_SK_SZ + MM_PK_SZ)];

#ifdef MM_KEM
    uint8_t ct[MM_CTU_SZ + (MM_N_MAX * MMKEM_CTI_SZ)];
    uint8_t kk[MM_N_MAX * MMKEM_K_SZ];
#else
    uint8_t ct[MM_CTU_SZ + (MM_N_MAX * MMPKE_CTI_SZ)];
    uint8_t mm[MM_N_MAX * MMPKE_M_SZ];
#endif

    printf( "%16s  %16s  m= %d  n= %d  du= %d  dv= %d  sig0= %f  sig1= %f\n",
            MM_PAR, "parameters",
            MM_M, MM_N, MM_DU, MMPKE_DV, MM_SIGMA0, MM_SIGMA1);

#ifdef  TESTVEC
    nn      =   5;
    rep     =   1;
    printf("=== chk %s N=%d\n", MM_PAR, nn);
#else
    int seed = time(NULL);

    printf( "%16s  %16s  time= %d\n",
            MM_PAR, "seed", seed);
    srandom(seed);
    seed = random();

    //  nn recipients loop
#ifdef  MM_SET_N
    nn  =   MM_SET_N; {                     //  set specific N
#else
    for (nn = 1; nn <= 1024; nn *= 2) {
#endif

    for (i = 0; i < 16; i++) {
        seed_a[i] = random();
    }
    for (i = 0; i < 32; i++) {
        seed_k[i] = random();
        seed_e[i] = random();
    }

    rep =   MM_REP_TOT / nn;

#endif

#ifdef MM_KEM
    nn_ct_sz = MM_CTU_SZ + nn * MMKEM_CTI_SZ;
#else
    nn_ct_sz = MM_CTU_SZ + nn * MMPKE_CTI_SZ;
#endif
    printf( "%16s  %16s  N= %4d  len= %9zu\n",
            MM_PAR, "ciphertext", nn, nn_ct_sz);

    //  --- mmSetup() ---

    dd  = get_sec();
    cc  = plat_get_cycle();
    for (iter = 0; iter < rep; iter++) {
        mm_setup(a_mat, seed_a);
    }
    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;

    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            MM_PAR, "mmSetup()", nn, cc, dd);

    //  --- mmKGen() ---

    dd  = get_sec();
    cc  = plat_get_cycle();
    for (iter = 0; iter < rep; iter++) {

        p = kdata;
        for (i = 0; i < nn; i++) {
            put64u_le(seed_k, i);
            sk[i] = p;
            p += MM_SK_SZ;
            pk[i] = p;
            mm_kgen(pk[i], sk[i], a_mat, seed_k);
            p += MM_PK_SZ;
#ifdef TESTVEC
            dbg_sum(sk[i], MM_SK_SZ, "sk");
            dbg_sum(pk[i], MM_PK_SZ, "pk");
#endif
        }
    }
    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            MM_PAR, "mmKGen()", nn, cc, dd);

#ifdef MM_PKE
    //  create random messages
    p   = mm;
    for (i = 0; i < nn; i++) {
        memcpy(buf, seed_e, 32);
        put64u_le(buf + 32, i);
        buf[40] = 'm';
        shake128(p, MMPKE_M_SZ, buf, 41);
        p += MMPKE_M_SZ;
    }
#endif


    dd  = get_sec();
    cc  = plat_get_cycle();
#ifdef MM_KEM
    //  --- mmEncap() ---
    for (iter = 0; iter < rep; iter++) {
        mm_encap(ct, kk, a_mat, (const uint8_t **) pk, seed_e, nn);
    }
#else
    //  --- mmEnc() ---
    for (iter = 0; iter < rep; iter++) {
        mm_enc(ct, a_mat, (const uint8_t **) pk, mm, seed_e, nn);
    }
#endif
    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;

    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
#ifdef MM_KEM
            MM_PAR, "mmEncap()", nn, cc, dd);
#else
            MM_PAR, "mmEnc()", nn, cc, dd);
#endif

#ifdef TESTVEC
    dbg_sum(ct, MM_CTU_SZ, "ct_u");
    dbg_sum(ct, nn_ct_sz, "ct");
#endif

    dd  = get_sec();
    cc  = plat_get_cycle();

    for (iter = 0; iter < rep; iter++) {
        for (i = 0; i < nn; i++) {
#ifdef MM_KEM
            //  --- mmDecap() ---
            uint8_t ki[MMKEM_K_SZ];
            uint8_t *cti = ct + MM_CTU_SZ + (i * MMKEM_CTI_SZ);

            mm_decap(ki, sk[i], ct, cti);

            //  report failures even in benchmarking (also, use result)
            if (memcmp(ki, kk + (i * MMKEM_K_SZ), MMKEM_K_SZ) != 0) {
                printf("[FAIL] decaps #%d\n", i);
            }
#ifdef TESTVEC
            printf("decaps #%d\n", i);
            dbg_sum(ki, MMKEM_K_SZ, "k_i");
            dbg_sum(cti, MMKEM_CTI_SZ, "ct_i");
#endif
#else
            //  --- mmDec() ---
            uint8_t mi[MMPKE_M_SZ];
            uint8_t *cti = ct + MM_CTU_SZ + (i * MMPKE_CTI_SZ);

            mm_dec(mi, sk[i], ct, cti);

            //  report failures even in benchmarking (also, use result)
            if (memcmp(mi, mm + (i * MMPKE_M_SZ), MMPKE_M_SZ) != 0) {
                printf("[FAIL] dec #%d\n", i);
            }
#ifdef TESTVEC
            printf("dec #%d\n", i);
            dbg_sum(mi, MMPKE_M_SZ, "m_i");
            dbg_sum(cti, MMPKE_CTI_SZ, "ct_i");
#endif
#endif
        }
    }
    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
#ifdef MM_KEM
            MM_PAR, "mmDecap()", nn, cc, dd);
#else
            MM_PAR, "mmDec()", nn, cc, dd);
#endif

#ifndef TESTVEC
    }   //  nn recipients loop
#endif

    return 0;
}

#endif
