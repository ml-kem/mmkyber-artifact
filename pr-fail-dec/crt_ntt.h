//  crt_ntt.h
//  === Header: CRT NTT Convolution

#ifndef _CRT_NTT_H_
#define _CRT_NTT_H_

#include <gmp.h>

//  initialize the ntt operations (call once)
size_t crt_ntt_init();

//  free static structures
void crt_ntt_free();

//  compute a convolution r += f * g. r needs to hold "bits"
size_t crt_ntt_convol(  mpz_t r[],
                        mpz_t f[], size_t f_len,
                        mpz_t g[], size_t g_len,
                        size_t bits);

#endif
