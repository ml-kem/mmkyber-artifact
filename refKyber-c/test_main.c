//  test_main.c
//  === benchmark script of plain C Kyber

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "plat_local.h"
#include "kem.h"
#include "indcpa.h"

#ifndef MM_N_MAX
#define MM_N_MAX 1024
#endif

#ifndef MM_REP_TOT
#define MM_REP_TOT 10240
#endif

static double get_sec()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    return ((double) tv.tv_sec) + 1E-6*((double) tv.tv_usec);
}

int main()
{
    const size_t    max_n = MM_N_MAX;
    const size_t    pk_sz = CRYPTO_PUBLICKEYBYTES;
    const size_t    sk_sz = CRYPTO_SECRETKEYBYTES;
    const size_t    ct_sz = CRYPTO_CIPHERTEXTBYTES;
    const size_t    ss_sz = CRYPTO_BYTES;
    const char      *algn = CRYPTO_ALGNAME;

    uint8_t seed_dz[64] = { 0 };
    uint8_t seed_m[32] = { 0 };
    uint8_t seed_r[32] = { 0 };

    uint8_t pk[max_n][pk_sz];
    uint8_t sk[max_n][sk_sz];
    uint8_t ct[max_n][ct_sz];
    uint8_t ss[max_n][ss_sz];
    uint8_t s2[ss_sz];

    uint64_t cc = 0;    //  cycle counts
    double dd = 0.0;    //  seconds

    int i, nn, rep, iter, seed;

    printf( "%16s  %16s  |pk|= %zu  |sk|= %zu  |ct|= %zu  |ss|= %zu\n",
            algn, "parameters", pk_sz, sk_sz, ct_sz, ss_sz);

    seed = time(NULL);
    printf( "%16s  %16s  time= %d\n",
            algn, "seed", seed);
    srandom(seed);
    seed    = random();

#ifdef  MM_SET_N
    nn  =   MM_SET_N; {                     //  set specific N
#else
    for (nn = 1; nn <= 1024; nn *= 2) {
#endif

    printf( "%16s  %16s  N= %4d  len= %9d\n",
            algn, "ciphertext", nn, (int) ct_sz * nn);

    //  this is mostly redudndant
    for (i = 0; i < 32; i++) {
        seed_dz[i]      = random();     //  keygen d
        seed_dz[i + 32] = random();     //  keygen z
        seed_m[i]       = random();     //  encaps m
        seed_r[i]       = random();     //  encrypt r
    }

    rep     = MM_REP_TOT / nn;

    //  --- key gen speed
    dd  = get_sec();
    cc  = plat_get_cycle();

    for (iter = 0; iter < rep; iter++) {
        for (i = 0; i < nn; i++) {
            put64u_le(seed_dz, seed++);         //  change d seed
            put64u_le(seed_dz + 32, seed++);    //  change z (not so relevant)

            //  equivalent to Algorithm 16, ML-KEM.KeyGen_internal(d, z)
            crypto_kem_keypair_derand(pk[i], sk[i], seed_dz);
        }
    }
    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            algn, "ML-KEM.KeyGen()", nn, cc, dd);

    //  --- encaps speed

    dd  = get_sec();
    cc  = plat_get_cycle();

    for (iter = 0; iter < rep; iter++) {
        for (i = 0; i < nn; i++) {
            put64u_le(seed_m, seed++);          //  change m seed

            //  equivalent to Algorithm 17, ML-KEM.Encaps_internal(ek, m)
            crypto_kem_enc_derand(ct[i], ss[i], pk[i], seed_m);
        }
    }

    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            algn, "ML-KEM.Encaps()", nn, cc, dd);

    //  --- decaps speed

    dd  = get_sec();
    cc  = plat_get_cycle();

    for (iter = 0; iter < rep; iter++) {
        for (i = 0; i < nn; i++) {

            //  equivalent to Algorithm 21, ML-KEM.Decaps(dk, c)
            crypto_kem_dec(s2, ct[i], sk[i]);

            if (memcmp(s2, ss[i], ss_sz) != 0) {
                printf("[FAIL] decaps #%d\n", i);
            }
        }
    }

    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            algn, "ML-KEM.Decaps()", nn, cc, dd);

    //  --- encrypt speed

    dd  = get_sec();
    cc  = plat_get_cycle();

    for (iter = 0; iter < rep; iter++) {
        for (i = 0; i < nn; i++) {
            put64u_le(seed_r, seed++);          //  change r seed

            //  equivalent to Algorithm 14, K-PKE.Encrypt(ekPKE, m, r)
            indcpa_enc(ct[i], ss[i], pk[i], seed_r);
            //  (we use the previous shared secrets as message inputs here)
        }
    }

    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            algn, "K-PKE.Encrypt()", nn, cc, dd);


    //  --- decrypt speed

    dd  = get_sec();
    cc  = plat_get_cycle();

    for (iter = 0; iter < rep; iter++) {
        for (i = 0; i < nn; i++) {

            //  equivalent to Algorithm 15, K-PKE.Decrypt(dkPKE , c)
            indcpa_dec(s2, ct[i], sk[i]);

            if (memcmp(s2, ss[i], ss_sz) != 0) {
                printf("[FAIL] decrypt #%d\n", i);
            }
        }
    }

    cc  = (plat_get_cycle() - cc) / rep;
    dd  = (get_sec() - dd) / rep;
    printf( "%16s  %16s  N= %4d  cyc= %9lu  sec= %8.6f\n",
            algn, "K-PKE.Decrypt()", nn, cc, dd);

    }   //  nn recipients loop

    return 0;
}
