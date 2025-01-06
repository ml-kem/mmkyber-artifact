//  mm_sample.h
//  === Header: General Polynomial Ring Arithmetic.

#ifndef _MM_SAMPLE_H_
#define _MM_SAMPLE_H_

#include "plat_local.h"
#include "sha3_t.h"

//  Uniform sampler [0, q-1]
void poly_unif(int32_t *r, sha3_t *kec);

//  Sample bytes for a nu-distribution polynomial
void sample_nu(uint8_t *r, sha3_t *kec);

//  Decode a nu-distribution polynomial from bytes
void poly_nu(int32_t *r, const uint8_t *s);

//  Gaussian sampler. gw = Gaussian width, gw = sqrt(2*Pi)*sigma
void poly_gauss(int32_t *r, sha3_t *kec, double gw);

#endif
