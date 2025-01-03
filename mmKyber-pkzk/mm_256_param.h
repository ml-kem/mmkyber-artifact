// auto-generated by lnp-tbox.sage.
//
// protocol is statistically complete with correctness error >= 1 - 2^(-4)
// protocol is simulatable under MLWE(31,26,[-1,1])
// protocol is knowledge-sound with knowledge error <= 2^(-127.0) under MSIS(14,129,2^(33.35157))
//
// Ring
// degree d = 64
// modulus q = 562949953422733, log(q) ~ 49.0
// factors q = q1
//
// Compression
// D = 5
// gamma = 13932, log(gamma) ~ 13.766115
//
// Dimensions of secrets
// s1: m1 = 72
// m: l = 0
// s2: m2 = 57
//
// Size of secrets
// l2(s1) <= alpha = 68.0
// m unbounded
// s2 uniform in [-nu,nu] = [-1,1]
//
// Norm proofs
// binary: yes (dimension: 72)
// exact euclidean: no
// approximate infinity: yes (psi: 2983.8063, dimension: 36, bound: 2304.4999)
//
// Challenge space
// c uniform in [-omega,omega] = [-8,8], o(c)=c, sqrt(l1(o(c)*c)) <= eta = 140
//
// Standard deviations
// stdev1 = 101580.8, log(stdev1/1.55) = 16.0
// stdev2 = 6348.8, log(stdev2/1.55) = 12.0
// stdev3 = 6348.8, log(stdev3/1.55) = 12.0
// stdev4 = 203161.6, log(stdev4/1.55) = 17.0
//
// Repetition rate
// M1 = 3.5172203
// M2 = 2.4277063
// M3 = 1.01945
// M4 = 1.0219173
// total = 8.8956433
//
// Security
// MSIS dimension: 14
// MSIS root hermite factor: 1.0043999
// MLWE dimension: 31
// MLWE root hermite factor: 1.0043649
//
// Proof size
// ~ 28.4570312500000 KiB
//
// 50 bit moduli for degree 64: [1125899906840833, 1125899906839937, 1125899906837633]
// bit length of products: [49, 99, 149]
// inverses: [1, -162099428551732, 296975494591860]

#include "lazer.h"
static const limb_t _mm_256_param_q_limbs[] = {562949953422733UL};
static const int_t _mm_256_param_q = {{(limb_t *)_mm_256_param_q_limbs, 1, 0}};
static const limb_t _mm_256_param_qminus1_limbs[] = {562949953422732UL};
static const int_t _mm_256_param_qminus1 = {{(limb_t *)_mm_256_param_qminus1_limbs, 1, 0}};
static const limb_t _mm_256_param_m_limbs[] = {40406973401UL};
static const int_t _mm_256_param_m = {{(limb_t *)_mm_256_param_m_limbs, 1, 0}};
static const limb_t _mm_256_param_mby2_limbs[] = {40406973401/2UL};
static const int_t _mm_256_param_mby2 = {{(limb_t *)_mm_256_param_mby2_limbs, 1, 0}};
static const limb_t _mm_256_param_gamma_limbs[] = {13932UL};
static const int_t _mm_256_param_gamma = {{(limb_t *)_mm_256_param_gamma_limbs, 1, 0}};
static const limb_t _mm_256_param_gammaby2_limbs[] = {6966UL};
static const int_t _mm_256_param_gammaby2 = {{(limb_t *)_mm_256_param_gammaby2_limbs, 1, 0}};
static const limb_t _mm_256_param_pow2D_limbs[] = {32UL};
static const int_t _mm_256_param_pow2D = {{(limb_t *)_mm_256_param_pow2D_limbs, 1, 0}};
static const limb_t _mm_256_param_pow2Dby2_limbs[] = {16UL};
static const int_t _mm_256_param_pow2Dby2 = {{(limb_t *)_mm_256_param_pow2Dby2_limbs, 1, 0}};
static const limb_t _mm_256_param_Bsq_limbs[] = {668892785239UL, 0UL};
static const int_t _mm_256_param_Bsq = {{(limb_t *)_mm_256_param_Bsq_limbs, 2, 0}};
static const limb_t _mm_256_param_scM1_limbs[] = {7355360059028332997UL, 9541030083301590535UL, 3UL};
static const int_t _mm_256_param_scM1 = {{(limb_t *)_mm_256_param_scM1_limbs, 3, 0}};
static const limb_t _mm_256_param_scM2_limbs[] = {11208912776058102363UL, 7889787913383536571UL, 2UL};
static const int_t _mm_256_param_scM2 = {{(limb_t *)_mm_256_param_scM2_limbs, 3, 0}};
static const limb_t _mm_256_param_scM3_limbs[] = {16817174804370679746UL, 358788469537454048UL, 1UL};
static const int_t _mm_256_param_scM3 = {{(limb_t *)_mm_256_param_scM3_limbs, 3, 0}};
static const limb_t _mm_256_param_scM4_limbs[] = {7305702446238145672UL, 404302244700510563UL, 1UL};
static const int_t _mm_256_param_scM4 = {{(limb_t *)_mm_256_param_scM4_limbs, 3, 0}};
static const limb_t _mm_256_param_stdev1sq_limbs[] = {10318658929UL, 0UL};
static const int_t _mm_256_param_stdev1sq = {{(limb_t *)_mm_256_param_stdev1sq_limbs, 2, 0}};
static const limb_t _mm_256_param_stdev2sq_limbs[] = {40307261UL, 0UL};
static const int_t _mm_256_param_stdev2sq = {{(limb_t *)_mm_256_param_stdev2sq_limbs, 2, 0}};
static const limb_t _mm_256_param_stdev3sq_limbs[] = {40307261UL, 0UL};
static const int_t _mm_256_param_stdev3sq = {{(limb_t *)_mm_256_param_stdev3sq_limbs, 2, 0}};
static const limb_t _mm_256_param_stdev4sq_limbs[] = {41274635715UL, 0UL};
static const int_t _mm_256_param_stdev4sq = {{(limb_t *)_mm_256_param_stdev4sq_limbs, 2, 0}};
static const limb_t _mm_256_param_inv2_limbs[] = {281474976711366UL};
static const int_t _mm_256_param_inv2 = {{(limb_t *)_mm_256_param_inv2_limbs, 1, 1}};
static const limb_t _mm_256_param_inv4_limbs[] = {140737488355683UL};
static const int_t _mm_256_param_inv4 = {{(limb_t *)_mm_256_param_inv4_limbs, 1, 1}};
static const unsigned int _mm_256_param_n[0] = {};
static const limb_t _mm_256_param_Bz3sqr_limbs[] = {27753065054UL, 0UL};
static const int_t _mm_256_param_Bz3sqr = {{(limb_t *)_mm_256_param_Bz3sqr_limbs, 2, 0}};
static const limb_t _mm_256_param_Bz4_limbs[] = {3250585UL};
static const int_t _mm_256_param_Bz4 = {{(limb_t *)_mm_256_param_Bz4_limbs, 1, 0}};
static const limb_t _mm_256_param_Pmodq_limbs[] = {200649007881UL};
static const int_t _mm_256_param_Pmodq = {{(limb_t *)_mm_256_param_Pmodq_limbs, 1, 1}};
static const limb_t _mm_256_param_Ppmodq_0_limbs[] = {43308657UL};
static const int_t _mm_256_param_Ppmodq_0 = {{(limb_t *)_mm_256_param_Ppmodq_0_limbs, 1, 0}};
static const limb_t _mm_256_param_Ppmodq_1_limbs[] = {36290289UL};
static const int_t _mm_256_param_Ppmodq_1 = {{(limb_t *)_mm_256_param_Ppmodq_1_limbs, 1, 0}};
static const limb_t _mm_256_param_Ppmodq_2_limbs[] = {25615857UL};
static const int_t _mm_256_param_Ppmodq_2 = {{(limb_t *)_mm_256_param_Ppmodq_2_limbs, 1, 0}};
static const int_srcptr _mm_256_param_l2Bsq[] = {};
static const int_srcptr _mm_256_param_Ppmodq[] = {_mm_256_param_Ppmodq_0, _mm_256_param_Ppmodq_1, _mm_256_param_Ppmodq_2};
static const polyring_t _mm_256_param_ring = {{_mm_256_param_q, 64, 50, 6, moduli_d64, 3, _mm_256_param_Pmodq, _mm_256_param_Ppmodq, _mm_256_param_inv2}};
static const dcompress_params_t _mm_256_param_dcomp = {{ _mm_256_param_q, _mm_256_param_qminus1, _mm_256_param_m, _mm_256_param_mby2, _mm_256_param_gamma, _mm_256_param_gammaby2, _mm_256_param_pow2D, _mm_256_param_pow2Dby2, 5, 1, 36 }};
static const abdlop_params_t _mm_256_param_tbox = {{ _mm_256_param_ring, _mm_256_param_dcomp, 72, 57, 0, 12, 14, _mm_256_param_Bsq, 1, 8, 5, 140, 1, 16, _mm_256_param_scM1, _mm_256_param_stdev1sq, 2, 12, _mm_256_param_scM2, _mm_256_param_stdev2sq}};
static const abdlop_params_t _mm_256_param_quad_eval_ = {{ _mm_256_param_ring, _mm_256_param_dcomp, 72, 57, 9, 3, 14, _mm_256_param_Bsq, 1, 8, 5, 140, 1, 16, _mm_256_param_scM1, _mm_256_param_stdev1sq, 2, 12, _mm_256_param_scM2, _mm_256_param_stdev2sq}};
static const abdlop_params_t _mm_256_param_quad_many_ = {{ _mm_256_param_ring, _mm_256_param_dcomp, 72, 57, 11, 1, 14, _mm_256_param_Bsq, 1, 8, 5, 140, 1, 16, _mm_256_param_scM1, _mm_256_param_stdev1sq, 2, 12, _mm_256_param_scM2, _mm_256_param_stdev2sq}};
static const lnp_quad_eval_params_t _mm_256_param_quad_eval = {{ _mm_256_param_quad_eval_, _mm_256_param_quad_many_, 4}};
static const lnp_tbox_params_t _mm_256_param = {{ _mm_256_param_tbox, _mm_256_param_quad_eval, 72, _mm_256_param_n, 36, 0, 72, 2, 12, _mm_256_param_scM3, _mm_256_param_stdev3sq, 2, 17, _mm_256_param_scM4, _mm_256_param_stdev4sq, _mm_256_param_Bz3sqr, _mm_256_param_Bz4, &_mm_256_param_l2Bsq[0], _mm_256_param_inv4, 29140UL }};

static const unsigned int mm_256_param_Ps[72] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71};

static const limb_t mm_256_param_p_limbs[] = {33550337UL};
static const int_t mm_256_param_p = {{(limb_t *)mm_256_param_p_limbs, 1, 0}};
static const limb_t mm_256_param_pinv_limbs[] = {34187934462426UL};
static const int_t mm_256_param_pinv = {{(limb_t *)mm_256_param_pinv_limbs, 1, 1}};
static const unsigned int mm_256_param_s1_indices[18] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
static const lin_params_t mm_256_param = {{ _mm_256_param, 256, mm_256_param_p, mm_256_param_pinv, 4, mm_256_param_s1_indices, 18, NULL, 0,  mm_256_param_Ps, 72, NULL, NULL, NULL, NULL }};

