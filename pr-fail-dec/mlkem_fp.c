//  mlkem_fp.c
//  === ML-KEM Failure Probability estimate

#include <stdio.h>
#include <math.h>
#include "zp_dist.h"

//  derived from Theorem 1 in original Kyber paper

double mlkem_fp(    int64_t     prm_k,
                    int64_t     prm_eta1,
                    int64_t     prm_eta2,
                    int64_t     prm_d_u,
                    int64_t     prm_d_v,
                    size_t      eps)
{
    const int64_t   prm_q       = 3329;     //  13*256+1
    const int64_t   prm_n       = 256;      //  polynomial degree
    double var, tail = 0.0;

    printf("\n=== mlkem_fp( %ld, %ld, %ld, %ld, %ld, %zu )\n",
            prm_k, prm_eta1, prm_eta2, prm_d_u, prm_d_v, eps);

    //  temporary
    dist_t *ds = NULL, *d1 = NULL, *d2 = NULL, *d3 = NULL;

    dist_init(eps);

    ds  = dist_alloc(1, eps);   //  sum
    d1  = dist_alloc(1, eps);   //  temporary 1-3
    d2  = dist_alloc(1, eps);
    d3  = dist_alloc(1, eps);

    //  === <e, y>
    dist_cbd(   d1, prm_eta1, eps);
    dist_print( d1, "e");

    dist_cbd(   d2, prm_eta1, eps);
    dist_print( d2, "y");

    dist_mul(   d3, d1, d2);
    dist_print( d3, "e*y");

    dist_scmul( ds, d3, prm_k * prm_n);
    dist_print( ds, "<e,y>");


    //  === e2
    dist_cbd(   d1, prm_eta2, eps);
    dist_print( d1, "e2");

    dist_add(   ds, ds, d1);
    dist_print( ds, "SUM: <e,y> + e2");


    //  === c_v
    dist_round( d1, prm_d_v, prm_q, eps);
    dist_print( d1, "c_v");

    dist_add(   ds, ds, d1);
    dist_print( ds, "SUM: <e,y> + e2 + c_v");


    //  === <s, e_1>
    dist_cbd(   d1, prm_eta1, eps);
    dist_print( d1, "s");

    dist_cbd(   d2, prm_eta2, eps);
    dist_print( d2, "e_1");

    dist_mul(   d3, d1, d2);
    dist_print( d3, "s*e_1");

    dist_scmul( d1, d3, prm_k * prm_n);
    dist_print( d1, "<s,e_1>");

    dist_neg(   d2, d1  );
    dist_add(   ds, ds, d2);
    dist_print( ds, "SUM: <e,y> + e2 + c_v - <s,e_1>");


    //  === <s, c_u>
    dist_cbd(   d1, prm_eta1, eps);
    dist_print( d1, "s");

    dist_round( d2, prm_d_u, prm_q, eps);
    dist_print( d2, "c_u");

    dist_mul(   d3, d1, d2);
    dist_print( d3, "s*c_u");

    dist_scmul( d1, d3, prm_k * prm_n);
    dist_print( d1, "<s,c_u>");

    dist_neg(   d2, d1);
    dist_add(   ds, ds, d2);
    var = dist_print( ds, "SUM: <e,y> + e2 + c_v - <s,e_1> - <s,c_u>");
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

double mlkem_fp_summary(size_t eps)
{
    //  the probabilities are so small that we may approximate the product
    //  probability  1 - prod_i^n (1-p_i) as n*p_i

    double fp_128 = 256.0 * mlkem_fp(2, 3, 2, 10, 4, eps);
    double fp_192 = 256.0 * mlkem_fp(3, 2, 2, 10, 4, eps);
    double fp_256 = 256.0 * mlkem_fp(4, 2, 2, 11, 5, eps);

    printf("\n=== ML-KEM Summary\n");
    printf("ML-KEM-512:  %e  2^%g\n", fp_128, log2(fp_128));
    printf("ML-KEM-768:  %e  2^%g\n", fp_192, log2(fp_192));
    printf("ML-KEM-1024: %e  2^%g\n", fp_256, log2(fp_256));

    double fp = fp_128 < fp_192 ? fp_128 : fp_192;
    fp = fp < fp_256 ? fp : fp_256;

    return fp;
}
