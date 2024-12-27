//  mmkyber.h
//  === Header: mmKyber-KEM and mmKyber-PKE implemetation ("public API")

#ifndef _MMKYBER_H_
#define _MMKYBER_H_

#include <stdint.h>
#include <stddef.h>

//  mmKEM & mmPKE: mmSetup(1^lambda, N): Generate public parameter A from seed.
void mm_setup(  int32_t *a, const uint8_t seed_a[16]);

//  mmKEM & mmPKE: mmKGen(pp): Generate individual public key.
size_t mm_kgen( uint8_t *pk, uint8_t *sk,
                const int32_t *a_mat, const uint8_t seed_k[32]);

//  mmKEM: mmEncap(pp, (pk_i) for i in [N]): Encapsulate to N recipients.
size_t mm_encap(uint8_t *ct, uint8_t *kk,
                const int32_t *a_mat, const uint8_t *pk[],
                const uint8_t seed_e[32], size_t n);

//  mmKEM: mmDecap(pp, sk, ct): Decapsulate individual ciphertext (ctu,cti).
void mm_decap(  uint8_t *k, const uint8_t *sk,
                const uint8_t *ctu, const uint8_t *cti);

//  mmPKE: mmEnc(pp, (pk_i), (m_i) for i in [N]): Encrypt to N recipients.
size_t mm_enc(  uint8_t *ct, const int32_t *a_mat,
                const uint8_t *pk[], const uint8_t *mm,
                const uint8_t seed_e[32], size_t n);

//  mmPKE: mmDec(pp, sk, ct): Decrypt a message
void mm_dec(uint8_t *m, const uint8_t *sk,
            const uint8_t *ctu, const uint8_t *cti);

#endif
