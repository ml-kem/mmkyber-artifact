//  mmkyber_fp.c
//  === mmKyber Failure Probability estimate

#include <stdio.h>
#include <math.h>
#include "zp_dist.h"

//  mmKyber Theorem D.6. (Correctness)

double mmkyber_fp(  int64_t     prm_n,
                    int64_t     prm_m,
                    int64_t     prm_bar_nu,
                    int64_t     prm_d_u,
                    int64_t     prm_d_v,
                    const char  *prm_sigma0,
                    const char  *prm_sigma1,
                    bool        kem,
                    size_t      eps)
{
    const int64_t   prm_q       = 33550337;     //  2**25 - 2**12 + 1
    const int64_t   prm_d       = 256;          //  polynomial degree
    double var, tail = 0.0;

    printf("\n=== mmkyber_fp( %ld, %ld, %ld, %ld, %ld, %s, %s, %s, %zu )\n",
            prm_n, prm_m, prm_bar_nu, prm_d_u, prm_d_v,
            prm_sigma0, prm_sigma1, kem ? "KEM" : "PKE", eps);

    //  temporary
    dist_t *ds = NULL, *d1 = NULL, *d2 = NULL, *d3 = NULL;

    //  initialize the NTT code (no harm if done multiple times)
    dist_init(eps);

    ds  = dist_alloc(1, eps);   //  sum
    d1  = dist_alloc(1, eps);   //  temporary 1-3
    d2  = dist_alloc(1, eps);
    d3  = dist_alloc(1, eps);

    //  === <e,r >
    dist_unif(  d1, -((prm_bar_nu - 1) / 2), prm_bar_nu / 2, eps);
    dist_print( d1, "e");

    dist_gauss( d2, "0", prm_sigma0, true, eps);
    dist_print( d2, "r");

    dist_mul(   d3, d1, d2);
    dist_print( d3, "e*r");

    //  multiply by d * n
    dist_scmul( ds, d3, prm_d * prm_n);
    dist_print( ds, "SUM: <e,r>");

    //  === <s_i,e_u>  (same as above..)
    dist_unif(  d1, -((prm_bar_nu - 1) / 2), prm_bar_nu / 2, eps);
    dist_print( d1, "s_i");

    dist_gauss( d2, "0", prm_sigma0, true, eps);
    dist_print( d2, "e_u");

    dist_mul(   d3, d1, d2);
    dist_print( d3, "s_i*e_u");

    //  multiply by d * m
    dist_scmul( d1, d3, prm_d * prm_m);
    dist_neg(   d2, d1  );
    dist_print( d2, "-<s_i,e_u>");

    dist_add(   ds, ds, d2);
    dist_print( ds, "SUM: <e,r> - <s_i,e_u>");


    //  === <s_i, c_u>
    dist_unif(  d1, -((prm_bar_nu - 1) / 2), prm_bar_nu / 2, eps);
    dist_print( d1, "s_i");

    dist_round( d2, prm_d_u, prm_q, eps);
    dist_print( d2, "c_u");

    dist_mul(   d3, d1, d2);
    dist_print( d3, "s_i*c_u");

    //  multiply by d * m
    dist_scmul( d1, d3, prm_d * prm_m);
    dist_print( d1, "<s_i,c_u>");

    dist_add(   ds, ds, d1);
    dist_print( ds, "SUM: <e,r> - <s_i,e_u> + <s_i,c_u>");

    //  (widest distributions are added last)
    //  === y_i
    dist_gauss( d1, "0", prm_sigma1, true, eps);
    dist_print( d1, "y_i");

    dist_add(   ds, ds, d1);
    dist_print( ds, "SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>");

    if (kem) {
        //  mmKyberKEM -- Theorem 4.2

        //  === scale by 2
        dist_spread(ds, ds, 2);
        dist_print( ds, "SUM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>)");

        //  === - c_v - e_i
        dist_dbl_r( d1, prm_d_u, prm_q, eps);
        dist_neg(   d2, d1  );
        dist_print( d2, "- c_v - e_i");

        dist_add(   ds, ds, d2);
        var = dist_print( ds,
            "KEM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>) - c_v - e_i");

    } else {

        //  mmKyberPKE -- Theorem 4.6

        //  === -c_v
        dist_round( d1, prm_d_v, prm_q, eps);
        dist_neg(   d2, d1  );
        dist_print( d2, "-c_v");

        dist_add(   ds, ds, d2);
        var = dist_print( ds,
            "PKE: <e,r> + y_i - <s_i,e_u> - c_v + <s_i,c_u>");
    }
    printf("Total sigma/q= %f\n", sqrt(var)/((double) prm_q));

    //  compute tail
    tail = dist_tail(ds, -prm_q/4, prm_q/4);

    printf("tail per coefficient = %e  2^%g\n", tail, log2(tail));

    dist_clear(ds);
    dist_clear(d1);
    dist_clear(d2);
    dist_clear(d3);

    return tail;
}

double mmkyber_fp_summary(size_t eps)
{
    //  the probabilities are so small that we may approximate the product
    //  probability  1 - prod_i^dn (1-p_i) as dn*p_i
    double fp_128, fp_192, fp_256;
    const int64_t n = 1024;
    const double dn = (double) (256 * n);   //  secret bits d * recipients N

    //  mmKyberKEM
    fp_128 = dn * mmkyber_fp(4, 4, 3, 10, 0, "15.90", "368459.34", true, eps);
    fp_192 = dn * mmkyber_fp(7, 7, 2, 11, 0, "15.90", "488797.36", true, eps);
    fp_256 = dn * mmkyber_fp(9, 9, 2, 11, 0, "15.90", "554941.07", true, eps);

    printf("\n=== mmKyberKEM Summary N= %ld\n", n);
    printf("mmKyberKEM-128: %e  2^%g\n", fp_128, log2(fp_128));
    printf("mmKyberKEM-192: %e  2^%g\n", fp_192, log2(fp_192));
    printf("mmKyberKEM-256: %e  2^%g\n", fp_256, log2(fp_256));


    //  mmKyberPKE
    fp_128 = dn * mmkyber_fp(4, 4, 3, 10, 2, "15.90", "368459.34", false, eps);
    fp_192 = dn * mmkyber_fp(7, 7, 2, 11, 2, "15.90", "488797.36", false, eps);
    fp_256 = dn * mmkyber_fp(9, 9, 2, 11, 2, "15.90", "554941.07", false, eps);

    printf("\n=== mmKyberPKE Summary N= %ld\n", n);
    printf("mmKyberPKE-128: %e  2^%g\n", fp_128, log2(fp_128));
    printf("mmKyberPKE-192: %e  2^%g\n", fp_192, log2(fp_192));
    printf("mmKyberPKE-256: %e  2^%g\n", fp_256, log2(fp_256));

    double fp = fp_128 < fp_192 ? fp_128 : fp_192;
    fp = fp < fp_256 ? fp : fp_256;

    return fp;
}

