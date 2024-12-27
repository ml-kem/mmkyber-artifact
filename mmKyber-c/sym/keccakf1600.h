//  keccakf1600.h
//  === Raw Keccak f-1600 interface.

#ifndef _KECCAKF1600_H_
#define _KECCAKF1600_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "plat_local.h"

//  == low-level interface, keccakf1600.c or native

//  FIPS 202 Keccak f1600 permutation, 24 rounds
void keccak_f1600(uint64_t state[25]);

//  clear the state
static inline void keccak_clear(uint64_t state[25])
{
    int i;

    for (i = 0; i < 25; i++) {
        state[i] = 0;
    }
}

//  extract "rate" bytes from state
static inline
void keccak_extract(uint64_t* state, uint8_t *data, int rate)
{
    int i;

    for (i = 0; i < rate / 8; i++) {
        put64u_le(data + 8 * i, state[i]);
    }
}

//  absorb "rate" bytes via xor into the state
static inline
void keccak_xorbytes(uint64_t* state, const uint8_t *data, int rate)
{
    int i;

    for (i = 0; i < rate / 8; i++) {
        state[i] ^= get64u_le(data + 8 * i);
    }
}


#ifdef __cplusplus
}
#endif

#endif  //  _KECCAKF1600_H_
