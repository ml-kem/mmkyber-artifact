//  mm_param.h
//  === Header: mmKyber-KEM and mmKyber-PKE parameter sets

#ifndef _MM_PARAM_H_
#define _MM_PARAM_H_

//  which one

#if defined(MM_KEM)
#define MM_NAME "mmKyber-KEM"
#elif defined(MM_PKE)
#define MM_NAME "mmKyber-PKE"
#else
#error "Neither MM_KEM or MM_PKE_defined."
#endif

#if defined(MM_128)

//  mmKyber-KEM & mmKyber-PKE-128
#define MM_PAR      (MM_NAME "-128")
#define MM_M        4
#define MM_N        4
#define MM_DU       14
#define MMPKE_DV    2
#define MM_NU_BAR   3
#define MM_SIGMA0   15.90
#define MM_SIGMA1   368459.34
#elif defined(MM_192)

//  mmKyber-KEM & mmKyber-PKE-192
#define MM_PAR      (MM_NAME "-192")
#define MM_M        7
#define MM_N        7
#define MM_DU       15
#define MMPKE_DV    2
#define MM_NU_BAR   2
#define MM_SIGMA0   15.90
#define MM_SIGMA1   488797.36

#elif defined(MM_256)

//  mmKyber-KEM & mmKyber-PKE-256
#define MM_PAR      (MM_NAME "-256")
#define MM_M        9
#define MM_N        9
#define MM_DU       16
#define MMPKE_DV    2
#define MM_NU_BAR   2
#define MM_SIGMA0   15.90
#define MM_SIGMA1   554941.07

#else
#error "Parameter set not set: One of MM_128, MM_192, MM_256."
#endif

//  small uniform samplers
#if (MM_NU_BAR == 2)
//  binary secret
#define MM_NU_SZ    (MM_D / 8)
#elif (MM_NU_BAR == 3)
//  ternary secret
#define MM_NU_SZ    ((MM_D + 4) / 5)
#else
#error "Parameter MM_NU_BAR not set.."
#endif

//  shared parameters
#define MM_D        256         //  base polynomial x^d + 1
#define MM_Q        33550337l   //  2**25 - 2**12 + 1
#define MM_LOGQ     25

//  derived parameters
#define MMKEM_K_SZ      32
#define MMPKE_M_SZ      32
#define MMKEM_CTI_SZ    (MM_D / 8)
#define MMPKE_CTI_SZ    ((MMPKE_DV * MM_D) / 8)

#define MM_CTU_SZ       ((MM_M * MM_DU * MM_D) / 8)
#define MM_PK_SZ        ((MM_N * MM_LOGQ * MM_D) / 8)
#define MM_SK_SZ        (MM_NU_SZ * MM_M)

#endif
