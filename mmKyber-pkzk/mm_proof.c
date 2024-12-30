//  mm_proof.c
//  === ZK PoK for mmKyber keys

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

//  from mmKyber-c directory
#include "mmkyber.h"
#include "mm_ring.h"
#include "mm_param.h"
#include "mm_sample.h"
#include "mm_serial.h"

//  from lazer
#include "lazer.h"

//  LaZer parameter sets

#if defined(MM_128)
#include "mm_128_param.h"
#define LIN_PARAM mm_128_param
#elif defined(MM_192)
#include "mm_192_param.h"
#define LIN_PARAM mm_192_param
#elif defined(MM_256)
#include "mm_256_param.h"
#define LIN_PARAM mm_256_param
#endif


void debug_key(int64_t *mm_a, int64_t *mm_se, int64_t *mm_t, int64_t q);

int64_t check_key(  const int64_t *mm_a, const int64_t *mm_se,
                    const int64_t *mm_t, int64_t q);

//  checksum

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


//  for indexing the A matrix

#define MM_A_IDX(i,j) (((i) * MM_N + (j)) * MM_D)

//  extract A matrix from NTT domain representation

void mm_extract_a(int64_t *p_a, const int32_t *a_mat)
{
    int i, j, k;
    int32_t a[MM_D];

    //  transpose and serialize A
    for (i = 0; i < MM_N; i++) {
        for (j = 0; j < MM_M; j++) {

            polyr_copy(a, &a_mat[MM_A_IDX(j, i)]);
            polyr_scale(a, 1, a);       //  Montgomery scaling
            polyr_intt(a);              //  inverse transform

            for (k = 0; k < MM_D; k++) {
                *p_a++ = (int64_t) a[k];
            }
        }
    }
}

//  extract (s | e) vector and t from serialized keypair (pk, sk)

void mm_extract_keys(   int64_t *p_se, int64_t *p_t,
                        const uint8_t *pk, const uint8_t *sk,
                        const int32_t *a_mat)
{
    int i, j;
    int32_t e[MM_D], t[MM_D];
    int32_t s[MM_M][MM_D];

    //  expand private key
    for (i = 0; i < MM_M; i++) {
        poly_nu(s[i], sk);
        sk += MM_NU_SZ;
        for (j = 0; j < MM_D; j++) {
            *p_se++ = (int64_t) s[i][j];
        }
        polyr_fntt(s[i]);
    }

    for (i = 0; i < MM_N; i++) {

        //  compute A*s
        polyr_zero(e);
        for (j = 0; j < MM_M; j++) {
            polyr_ntt_mul_add(e, &a_mat[MM_A_IDX(j, i)], s[j]);
        }
        polyr_intt(e);

        //  decode public key (in NTT domain)
        pk += poly_deserial(t, pk, MM_LOGQ);
        polyr_scale(t, 1, t);       //  Montgomery scaling
        polyr_intt(t);

        //  e = t - As
        polyr_sub(e, t, e);

        for (j = 0; j < MM_D; j++) {
            *p_se++ = (int64_t) e[j];
            *p_t++ = (int64_t) t[j];
        }
    }
}

//  gneric prover functions

size_t  prover( uint8_t *proof, polyvec_t s, polymat_t A, polyvec_t t,
                const uint8_t pp[32])
{
    size_t len = 0;

    lin_prover_state_t prover;

    lin_prover_init (prover, pp, LIN_PARAM);
    lin_prover_set_statement (prover, A, t);
    lin_prover_set_witness (prover, s);
    lin_prover_prove (prover, proof, &len, NULL);
    lin_prover_clear (prover);

    return len;
}

int verifier (  const uint8_t *proof, polymat_t A, polyvec_t t,
                const uint8_t pp[32])
{
    lin_verifier_state_t verifier;

    lin_verifier_init (verifier, pp, LIN_PARAM);
    lin_verifier_set_statement (verifier, A, t);

    int accept = lin_verifier_verify (verifier, proof, NULL);
    lin_verifier_clear (verifier);

    return accept;
}

//  helper: negate a (sub) matrix

void my_polymat_neg_self (polymat_t m)
{
    unsigned int i, j;
    poly_ptr ri;

    _MAT_FOREACH_ELEM (m, i, j)
    {
        ri = polymat_get_elem (m, i, j);
        poly_neg_self (ri);
    }
}

//  run standard key generation and extract vectors for LaZer

int mm_extract_tv(int64_t *p_a, int64_t *p_se, int64_t *p_t, int tv)
{
    //  A matrix
    int32_t a_mat[MM_M * MM_N * MM_D];              //  shared public paramter

    //  seeds
    uint8_t seed_a[] = "0123456789abcdef";
    uint8_t seed_k[] = "000102030405060708090a0b0c0d0e0f";

    //  serialized keypair
    uint8_t pk[MM_PK_SZ], sk[MM_SK_SZ] = { 0 };   //  public key, secret key

    //  set up A matrix
    mm_setup(a_mat, seed_a);

    //  extract for lazer
    mm_extract_a(p_a, a_mat);

    //  generate a keypair  --  standard key generation function
    put64u_le(seed_k, tv);  //  set test vector number
    mm_kgen( pk, sk, a_mat, seed_k );

    //  serialized public keys
    dbg_sum(pk, MM_PK_SZ, "pk");
    dbg_sum(sk, MM_SK_SZ, "sk");

    //  extract for lazer
    mm_extract_keys( p_se, p_t, pk, sk, a_mat);

    return 0;
}



#if (MM_NU_BAR == 3)

//  convert ternary (s | e) two into binary vectors

void split_tri2bin(int64_t *sese,   const int64_t *se)
{
    int i;
    int64_t *spos = sese;
    int64_t *sneg = spos + (2 * MM_N * MM_D);   //  negative half

    //  convert (s | e) into two binary vectors
    for (i = 0; i < 2 * MM_N * MM_D; i++) {
        switch (se[i]) {
            case 1:         //  use A | I
                spos[i] = 1;
                sneg[i] = 0;
                break;

            case -1:        //  use -A | -I
                spos[i] = 0;
                sneg[i] = 1;
                break;

            default:        //  zero
                spos[i] = 0;
                sneg[i] = 0;
                break;
        }
    }
}

#endif

//  for benchmarking..

static double get_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return ((double) tv.tv_sec) + 1E-6*((double) tv.tv_usec);
}

int main (void)
{
    int i, j;
    uint8_t pp[32] = { 0 };                         //  public randomness
    uint8_t proof[100000] = { 0 };                  //  puffer for the proof
    size_t total_sz, proof_sz = 0;
    double tim, proof_tim = 0.0, verify_tim = 0.0;

    //  expanded keys + error
    int64_t mm_a[MM_M * MM_N * MM_D] = { 0 };       //  shared public matrix A
    int64_t mm_se[(MM_M + MM_N) * MM_D] = { 0 };    //  secret key + error
    int64_t mm_t[MM_N * MM_D] = { 0 };              //  personal public t

    //  Note: As we have Id matrices, I'm just using the "n*n" for convenience.
    polymat_t A;        //  public matrix
    polyvec_t s;        //  secret
    polyvec_t t;        //  t=As

    INT_T (p, 1)
    POLYRING_T (Rp, p, 256);

#if (MM_NU_BAR == 2)

    //  binary secret (or L2 norm): [ A | I ]
    polymat_t Ap, Ip;                       //  submatrices
    polymat_alloc (A, Rp, MM_N, 2 * MM_N);
    polyvec_alloc (s, Rp, 2 * MM_N);        //  secret vector

    //  two submatrices:         | column |
    polymat_get_submat (Ap, A, 0,        0, MM_N, MM_N, 1, 1);
    polymat_get_submat (Ip, A, 0,     MM_N, MM_N, MM_N, 1, 1);

#elif (MM_NU_BAR == 3)

    //  ternary -> two binary secrets: [ A | I | -A | -I ]
    int64_t sese[4 * MM_N * MM_D] = { 0 };  //  conversion to 2 x binary
    polymat_t Ap, An, Ip, In;               //  pos. and neg. submatrices
    polymat_alloc (A, Rp, MM_N, 4 * MM_N);
    polyvec_alloc (s, Rp, 4 * MM_N);        //  secret vector

    //  four submatrices:        | column |
    polymat_get_submat (Ap, A, 0,        0, MM_N, MM_N, 1, 1);
    polymat_get_submat (Ip, A, 0,     MM_N, MM_N, MM_N, 1, 1);
    polymat_get_submat (An, A, 0, 2 * MM_N, MM_N, MM_N, 1, 1);
    polymat_get_submat (In, A, 0, 3 * MM_N, MM_N, MM_N, 1, 1);

#else
#error  "MM_NU_BAR ?"
#endif

    polyvec_alloc (t, Rp, MM_N);    //  public key

    //  some infos
    printf("%s: (m=%d, n=%d, du=%d, dv=%d, sig0=%f, sig1=%f)\n",
            MM_PAR, MM_M, MM_N, MM_DU, MMPKE_DV,
            MM_SIGMA0, MM_SIGMA1);
    printf("PoK of s: As-t=0, s in {-1,0,+1}\n");

    lazer_init ();

    total_sz = 0;
    proof_tim = 0.0;
    verify_tim = 0.0;

    for (i = 0; i < 10; i++) {

        printf("=== test vector %d\n", i);

        //  set the data from a test vector
        mm_extract_tv(mm_a, mm_se, mm_t, i);

#if (MM_NU_BAR == 2)

        //  [A | I ]
        polymat_set_i64 (Ap, mm_a);             //   A
        polymat_set_one (Ip);                   //   I
        polyvec_set_coeffvec_i64 (s, mm_se);    //   s = (s | e)

#elif (MM_NU_BAR == 3)

        //  [ A | I | -A | -I ]
        polymat_set_i64 (Ap, mm_a);             //   A  (for s pos. binary)
        polymat_set_one (Ip);                   //   I  (for e pos. binary)
        polymat_set_i64 (An, mm_a);             //  -A  (for s neg. binary)
        my_polymat_neg_self (An);
        polymat_set_one (In);                   //  -I  (for e neg. binary)
        my_polymat_neg_self (In);
        split_tri2bin(sese, mm_se);             //  convert ternary to bin
        polyvec_set_coeffvec_i64 (s, sese);     //  ( s+ | e+ | s- | e- )

#else
#error  "MM_NU_BAR ?"
#endif

        polyvec_set_coeffvec_i64 (t, mm_t);     //  set t
        polyvec_neg_self (t);                   //  negate: A*s = -t

        memset(proof, 0, sizeof(proof));        //  clear proof space

        //  "public randomness"
        for (j = 0; j < (int) sizeof(pp); j++) {
            pp[j] = random() & 0xFF;
        }

        printf ("create proof ... ");
        tim = get_sec();
        proof_sz = prover (proof, s, A, t, pp);
        tim = get_sec() - tim;
        proof_tim += tim;
        total_sz += proof_sz;

        printf("[OK] proof: %zu bytes. Time: %.6f sec.\n", proof_sz, tim);

        printf ("verify proof ... ");

        tim = get_sec();
        int accept = verifier (proof, A, t, pp);
        tim = get_sec() - tim;
        verify_tim += tim;

        printf ("%s Time: %.6f sec.\n", accept ? "[OK]" : "[FAILED]", tim);

    }

    printf("%s avg. zk proof size: %zu bytes.\n",
            MM_PAR, total_sz / i);
    printf("%s avg. zk proof gen time: %.6f sec.\n",
            MM_PAR, proof_tim / ((double) i));
    printf("%s avg. zk proof vfy time: %.6f sec.\n",
            MM_PAR, verify_tim / ((double) i));

    polymat_free (A);
    polyvec_free (s);
    polyvec_free (t);
    return 0;
}
