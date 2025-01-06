#   pr-fail-dec: Decryption/Decapsulation Failure Probability

This repo implements fast(er) probability distribution convolutions for
the purpose of estimating failure probabilities and other statistical
features in lattice-based cryptography.


##  What are we doing here?

Regular floating point arithmetic is generally not appropriate for calculating
failure probabilities (or you have to be very careful when using it), as the
failure probabilities should typically be under 2<sup>-128</sup>.
A 64-bit floating point numbers only provide 52 bits of precision; adding
floating point numbers of different magnitudes can become very imprecise --
and one has to do that when handling these distributions.

In this implementation, the distributions are represented with high-precision
fixed-point arithmetic; the distributions sum up to 2<sup>eps</sup>, where
eps is the precision parameter, typically 200 or above.
Convolutions (needed when adding distributions together) are implemented
with CRT-NTT, which gives that bottleneck operation a quasi-linear complexity.
Note, however, that the CRT "base" currently only has 32 primes
(see `crt_ntt.c`), being able to handle 1000-bit products, so you can't go
above eps > 500.

See `zp_dist.h` for available primitives, and `mlkem_fp.c` and `mmkyber_fp.c`
for examples of conducting tail mass calculations.


##  ML-KEM Failure Probability

Everyone seems to have their own numbers for Kyber decapsulation failure probability. Since they're all under 2<sup>-128</sup>, it's not so important in practice but it's weird that [FIPS 203](https://doi.org/10.6028/NIST.FIPS.203) has different numbers from the submitted Kyber specs, and the numbers given out by the referenced [official security scripts](https://github.com/pq-crystals/security-estimates).

The formulas used `mlkem_fp.c` in this repo have been adopted from Theorem 1 in the original [Kyber paper](https://eprint.iacr.org/2017/634), and they are essentially the same as in those in the official security scripts (the difference is in arithmetic precision.) Note that there have been changes to Kyber since Theorem 1 was published, such as the dropping of the "public key compression" parameter t; hence, c<sub>t</sub> is gone.

That security script from Léo Ducas and John Schanck uses 64-bit floating points (i.e., 53 bits of precision), so it is not surprising that it is a little off, even when they were careful when adding these convolution summands together. They do stuff like pruning the tails, which also contributes to the few bits of difference.

| Scheme      | This code            | Léo & John         | FIPS 203 Table 1   |
|-------------|----------------------|--------------------|--------------------|
| ML-KEM-512  | 2<sup>-141.588</sup> | 2<sup>-139.1</sup> | 2<sup>-138.8</sup> |
| ML-KEM-768  | 2<sup>-168.245</sup> | 2<sup>-165.2</sup> | 2<sup>-164.8</sup> |
| ML-KEM-1024 | 2<sup>-176.625</sup> | 2<sup>-175.2</sup> | 2<sup>-174.8</sup> |


##  mmKyberKEM and mmKyberPKE Decryption Failure Probability

The file `mmkyber_fp.c` implements the formulas of Theorem 4.2
(mmKyber-KEM correctness) and Theorem 4.6 (mmKyber-PKE correctness)
of the mmKyber paper.

Note that failure probability is for N=1024 recipients; all d=256 bits in all
of these N=1024 messages need to be decrypted correctly. As a close
approximation the the per-coefficient failure probability is scaled up by
factor 2<sup>18</sup>.

| Security   |   N  |      mmKyberKEM      |      mmKyberPKE      |
|------------|------|----------------------|----------------------|
| 128, Cat 1 | 1024 | 2<sup>-139.936</sup> | 2<sup>-148.856</sup> |
| 192, Cat 3 | 1024 | 2<sup>-204.656</sup> | 2<sup>-213.999</sup> |
| 256, Cat 5 | 1024 | 2<sup>-155.785</sup> | 2<sup>-164.723</sup> |


##  Building

The code was developed for a general Linux-like target. There is a
dependency on gmp and mpfr multiprecision libraries. Running
`sudo apt install libgmp-dev libmpfr-dev` takes care of that on Debian and
Ubuntu targets.

A basic `Makefile` is provided; this should generate an executable `xtest`,
which displays decryption failure estimates for ML-KEM, mmKyberKEM, and
mmKyberPKE.

```
=== mlkem_fp( 2, 3, 2, 10, 4, 256 )
e: avg= 0, std= 1.22474  [-3,3] (7)
y: avg= 0, std= 1.22474  [-3,3] (7)
e*y: avg= 0, std= 1.5  [-9,9] (19)
<e,y>: avg= 0, std= 33.9411  [-687,687] (2432)
e2: avg= 0, std= 1  [-2,2] (7)
SUM: <e,y> + e2: avg= 0, std= 33.9559  [-687,687] (2432)
c_v: avg= -0.0312406, std= 60.0638  [-104,104] (211)
SUM: <e,y> + e2 + c_v: avg= -0.0312406, std= 68.9976  [-781,781] (2432)
s: avg= 0, std= 1.22474  [-3,3] (211)
e_1: avg= 0, std= 1  [-2,2] (7)
s*e_1: avg= 0, std= 1.22474  [-6,6] (19)
<s,e_1>: avg= 0, std= 27.7128  [-556,556] (1688)
SUM: <e,y> + e2 + c_v - <s,e_1>: avg= -0.0312406, std= 74.355  [-947,947] (2432)
s: avg= 0, std= 1.22474  [-3,3] (1688)
c_u: avg= -0.000600781, std= 0.96125  [-2,2] (1113)
s*c_u: avg= 0, std= 1.17729  [-6,6] (19)
<s,c_u>: avg= 0, std= 26.6389  [-528,528] (1688)
SUM: <e,y> + e2 + c_v - <s,e_1> - <s,c_u>: avg= -0.0312406, std= 78.9829  [-1075,1075] (2432)
Total sigma/q= 0.023726
tail per coefficient = 9.322515e-46  2^-149.588

=== mlkem_fp( 3, 2, 2, 10, 4, 256 )
e: avg= 0, std= 1  [-2,2] (5)
y: avg= 0, std= 1  [-2,2] (5)
e*y: avg= 0, std= 1  [-4,4] (9)
<e,y>: avg= 0, std= 27.7128  [-541,541] (1600)
e2: avg= 0, std= 1  [-2,2] (5)
SUM: <e,y> + e2: avg= 0, std= 27.7308  [-541,541] (1600)
c_v: avg= -0.0312406, std= 60.0638  [-104,104] (211)
SUM: <e,y> + e2 + c_v: avg= -0.0312406, std= 66.1564  [-637,637] (1600)
s: avg= 0, std= 1  [-2,2] (211)
e_1: avg= 0, std= 1  [-2,2] (5)
s*e_1: avg= 0, std= 1  [-4,4] (9)
<s,e_1>: avg= 0, std= 27.7128  [-541,541] (1688)
SUM: <e,y> + e2 + c_v - <s,e_1>: avg= -0.0312406, std= 71.7263  [-842,842] (3200)
s: avg= 0, std= 1  [-2,2] (1688)
c_u: avg= -0.000600781, std= 0.96125  [-2,2] (1083)
s*c_u: avg= 0, std= 0.96125  [-4,4] (9)
<s,c_u>: avg= 0, std= 26.6389  [-515,515] (1688)
SUM: <e,y> + e2 + c_v - <s,e_1> - <s,c_u>: avg= -0.0312406, std= 76.5134  [-988,988] (3200)
Total sigma/q= 0.022984
tail per coefficient = 8.808864e-54  2^-176.245

=== mlkem_fp( 4, 2, 2, 11, 5, 256 )
e: avg= 0, std= 1  [-2,2] (5)
y: avg= 0, std= 1  [-2,2] (5)
e*y: avg= 0, std= 1  [-4,4] (9)
<e,y>: avg= 0, std= 32  [-618,618] (2304)
e2: avg= 0, std= 1  [-2,2] (5)
SUM: <e,y> + e2: avg= 0, std= 32.0156  [-618,618] (2304)
c_v: avg= -0.0156203, std= 30.034  [-52,52] (107)
SUM: <e,y> + e2 + c_v: avg= -0.0156203, std= 43.8981  [-663,662] (2304)
s: avg= 0, std= 1  [-2,2] (107)
e_1: avg= 0, std= 1  [-2,2] (5)
s*e_1: avg= 0, std= 1  [-4,4] (9)
<s,e_1>: avg= 0, std= 32  [-618,618] (1712)
SUM: <e,y> + e2 + c_v - <s,e_1>: avg= -0.0156203, std= 54.3235  [-901,901] (2304)
s: avg= 0, std= 1  [-2,2] (1712)
c_u: avg= -0.000300391, std= 0.620323  [-1,1] (1237)
s*c_u: avg= 0, std= 0.620323  [-2,2] (9)
<s,c_u>: avg= 0, std= 19.8503  [-384,384] (1712)
SUM: <e,y> + e2 + c_v - <s,e_1> - <s,c_u>: avg= -0.0156203, std= 57.8366  [-975,975] (2304)
Total sigma/q= 0.017374
tail per coefficient = 2.644329e-56  2^-184.625

=== ML-KEM Summary
ML-KEM-512:  2.386564e-43  2^-141.588
ML-KEM-768:  2.255069e-51  2^-168.245
ML-KEM-1024: 6.769482e-54  2^-176.625

=== mmkyber_fp( 4, 4, 3, 10, 0, 15.90, 368459.34, KEM, 256 )
e: avg= 0, std= 0.816497  [-1,1] (3)
r: avg= 0, std= 6.34318  [-118,118] (239)
e*r: avg= 0, std= 5.17919  [-118,118] (237)
SUM: <e,r>: avg= 0, std= 165.734  [-3125,3125] (7584)
s_i: avg= 0, std= 0.816497  [-1,1] (3)
e_u: avg= 0, std= 6.34318  [-118,118] (239)
s_i*e_u: avg= 0, std= 5.17919  [-118,118] (237)
-<s_i,e_u>: avg= 0, std= 165.734  [-3125,3125] (6251)
SUM: <e,r> - <s_i,e_u>: avg= 0, std= 234.383  [-4377,4377] (15168)
s_i: avg= 0, std= 0.816497  [-1,1] (7584)
c_u: avg= -0.000488281, std= 9458.15  [-16382,16382] (32767)
s_i*c_u: avg= 0, std= 7722.55  [-16382,16382] (32765)
<s_i,c_u>: avg= 0, std= 247122  [-4458058,4458058] (16775680)
SUM: <e,r> - <s_i,e_u> + <s_i,c_u>: avg= 0, std= 247122  [-4453019,4453019] (8906039)
y_i: avg= 0, std= 146994  [-2664412,2664412] (16775680)
SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>: avg= 0, std= 287535  [-5195038,5195038] (17812078)
SUM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>): avg= 0, std= 575070  [-10390076,10390076] (35624156)
- c_v - e_i: avg= 0.000976562, std= 18916.3  [-32764,32764] (65534)
KEM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>) - c_v - e_i: avg= 0.000976562, std= 575381  [-10366376,10366376] (35624156)
Total sigma/q= 0.017150
tail per coefficient = 2.861430e-48  2^-157.936

=== mmkyber_fp( 7, 7, 2, 11, 0, 15.90, 488797.36, KEM, 256 )
e: avg= 0.5, std= 0.5  [0,1] (2)
r: avg= 0, std= 6.34318  [-118,118] (239)
e*r: avg= 0, std= 4.48531  [-118,118] (237)
SUM: <e,r>: avg= 0, std= 189.872  [-3589,3589] (9104)
s_i: avg= 0.5, std= 0.5  [0,1] (2)
e_u: avg= 0, std= 6.34318  [-118,118] (239)
s_i*e_u: avg= 0, std= 4.48531  [-118,118] (237)
-<s_i,e_u>: avg= 0, std= 189.872  [-3589,3589] (7179)
SUM: <e,r> - <s_i,e_u>: avg= 0, std= 268.52  [-5020,5020] (18208)
s_i: avg= 0.5, std= 0.5  [0,1] (9104)
c_u: avg= -0.000244141, std= 4729.08  [-8191,8191] (16385)
s_i*c_u: avg= -0.00012207, std= 3343.96  [-8191,8191] (16383)
<s_i,c_u>: avg= -0.21875, std= 141557  [-2579359,2579359] (7339200)
SUM: <e,r> - <s_i,e_u> + <s_i,c_u>: avg= -0.21875, std= 141557  [-2576621,2576621] (5153243)
y_i: avg= 0, std= 195002  [-3534604,3534604] (7339200)
SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>: avg= -0.21875, std= 240965  [-4368113,4368112] (10306486)
SUM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>): avg= -0.4375, std= 481930  [-8736226,8736224] (20612972)
- c_v - e_i: avg= 0.000488281, std= 9458.15  [-16382,16382] (32770)
KEM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>) - c_v - e_i: avg= -0.437012, std= 482023  [-8713416,8713416] (20612972)
Total sigma/q= 0.014367
tail per coefficient = 9.412842e-68  2^-222.656

=== mmkyber_fp( 9, 9, 2, 11, 0, 15.90, 554941.07, KEM, 256 )
e: avg= 0.5, std= 0.5  [0,1] (2)
r: avg= 0, std= 6.34318  [-118,118] (239)
e*r: avg= 0, std= 4.48531  [-118,118] (237)
SUM: <e,r>: avg= 0, std= 215.295  [-4050,4050] (15168)
s_i: avg= 0.5, std= 0.5  [0,1] (2)
e_u: avg= 0, std= 6.34318  [-118,118] (239)
s_i*e_u: avg= 0, std= 4.48531  [-118,118] (237)
-<s_i,e_u>: avg= 0, std= 215.295  [-4050,4050] (8101)
SUM: <e,r> - <s_i,e_u>: avg= 0, std= 304.473  [-5676,5676] (15168)
s_i: avg= 0.5, std= 0.5  [0,1] (15168)
c_u: avg= -0.000244141, std= 4729.08  [-8191,8191] (16385)
s_i*c_u: avg= -0.00012207, std= 3343.96  [-8191,8191] (16383)
<s_i,c_u>: avg= -0.28125, std= 160510  [-2921171,2921170] (7766016)
SUM: <e,r> - <s_i,e_u> + <s_i,c_u>: avg= -0.28125, std= 160510  [-2918072,2918072] (5836145)
y_i: avg= 0, std= 221389  [-4004429,4004429] (15532032)
SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>: avg= -0.28125, std= 273454  [-4954579,4954578] (11672290)
SUM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>): avg= -0.5625, std= 546907  [-9909158,9909156] (23344580)
- c_v - e_i: avg= 0.000488281, std= 9458.15  [-16382,16382] (32770)
KEM: 2*(<e,r> + y_i - <s_i,e_u> + <s_i,c_u>) - c_v - e_i: avg= -0.562012, std= 546989  [-9882536,9882534] (23344580)
Total sigma/q= 0.016304
tail per coefficient = 4.848582e-53  2^-173.785

=== mmKyberKEM Summary N= 1024
mmKyberKEM-128: 7.501068e-43  2^-139.936
mmKyberKEM-192: 2.467520e-62  2^-204.656
mmKyberKEM-256: 1.271027e-47  2^-155.785

=== mmkyber_fp( 4, 4, 3, 10, 2, 15.90, 368459.34, PKE, 256 )
e: avg= 0, std= 0.816497  [-1,1] (3)
r: avg= 0, std= 6.34318  [-118,118] (239)
e*r: avg= 0, std= 5.17919  [-118,118] (237)
SUM: <e,r>: avg= 0, std= 165.734  [-3125,3125] (7584)
s_i: avg= 0, std= 0.816497  [-1,1] (3)
e_u: avg= 0, std= 6.34318  [-118,118] (239)
s_i*e_u: avg= 0, std= 5.17919  [-118,118] (237)
-<s_i,e_u>: avg= 0, std= 165.734  [-3125,3125] (6251)
SUM: <e,r> - <s_i,e_u>: avg= 0, std= 234.383  [-4377,4377] (15168)
s_i: avg= 0, std= 0.816497  [-1,1] (7584)
c_u: avg= -0.000488281, std= 9458.15  [-16382,16382] (32767)
s_i*c_u: avg= 0, std= 7722.55  [-16382,16382] (32765)
<s_i,c_u>: avg= 0, std= 247122  [-4458058,4458058] (16775680)
SUM: <e,r> - <s_i,e_u> + <s_i,c_u>: avg= 0, std= 247122  [-4453019,4453019] (8906039)
y_i: avg= 0, std= 146994  [-2664412,2664412] (16775680)
SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>: avg= 0, std= 287535  [-5195038,5195038] (17812078)
-c_v: avg= 0.125, std= 2.42129e+06  [-4193792,4193792] (8387585)
PKE: <e,r> + y_i - <s_i,e_u> - c_v + <s_i,c_u>: avg= 0.125, std= 2.4383e+06  [-9289036,9289036] (35624156)
Total sigma/q= 0.072676
tail per coefficient = 5.908565e-51  2^-166.856

=== mmkyber_fp( 7, 7, 2, 11, 2, 15.90, 488797.36, PKE, 256 )
e: avg= 0.5, std= 0.5  [0,1] (2)
r: avg= 0, std= 6.34318  [-118,118] (239)
e*r: avg= 0, std= 4.48531  [-118,118] (237)
SUM: <e,r>: avg= 0, std= 189.872  [-3589,3589] (9104)
s_i: avg= 0.5, std= 0.5  [0,1] (2)
e_u: avg= 0, std= 6.34318  [-118,118] (239)
s_i*e_u: avg= 0, std= 4.48531  [-118,118] (237)
-<s_i,e_u>: avg= 0, std= 189.872  [-3589,3589] (7179)
SUM: <e,r> - <s_i,e_u>: avg= 0, std= 268.52  [-5020,5020] (18208)
s_i: avg= 0.5, std= 0.5  [0,1] (9104)
c_u: avg= -0.000244141, std= 4729.08  [-8191,8191] (16385)
s_i*c_u: avg= -0.00012207, std= 3343.96  [-8191,8191] (16383)
<s_i,c_u>: avg= -0.21875, std= 141557  [-2579359,2579359] (7339200)
SUM: <e,r> - <s_i,e_u> + <s_i,c_u>: avg= -0.21875, std= 141557  [-2576621,2576621] (5153243)
y_i: avg= 0, std= 195002  [-3534604,3534604] (7339200)
SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>: avg= -0.21875, std= 240965  [-4368113,4368112] (10306486)
-c_v: avg= 0.125, std= 2.42129e+06  [-4193792,4193792] (8387585)
PKE: <e,r> + y_i - <s_i,e_u> - c_v + <s_i,c_u>: avg= -0.09375, std= 2.43325e+06  [-8475335,8475335] (20612972)
Total sigma/q= 0.072525
tail per coefficient = 1.450263e-70  2^-231.999

=== mmkyber_fp( 9, 9, 2, 11, 2, 15.90, 554941.07, PKE, 256 )
e: avg= 0.5, std= 0.5  [0,1] (2)
r: avg= 0, std= 6.34318  [-118,118] (239)
e*r: avg= 0, std= 4.48531  [-118,118] (237)
SUM: <e,r>: avg= 0, std= 215.295  [-4050,4050] (15168)
s_i: avg= 0.5, std= 0.5  [0,1] (2)
e_u: avg= 0, std= 6.34318  [-118,118] (239)
s_i*e_u: avg= 0, std= 4.48531  [-118,118] (237)
-<s_i,e_u>: avg= 0, std= 215.295  [-4050,4050] (8101)
SUM: <e,r> - <s_i,e_u>: avg= 0, std= 304.473  [-5676,5676] (15168)
s_i: avg= 0.5, std= 0.5  [0,1] (15168)
c_u: avg= -0.000244141, std= 4729.08  [-8191,8191] (16385)
s_i*c_u: avg= -0.00012207, std= 3343.96  [-8191,8191] (16383)
<s_i,c_u>: avg= -0.28125, std= 160510  [-2921171,2921170] (7766016)
SUM: <e,r> - <s_i,e_u> + <s_i,c_u>: avg= -0.28125, std= 160510  [-2918072,2918072] (5836145)
y_i: avg= 0, std= 221389  [-4004429,4004429] (15532032)
SUM: <e,r> + y_i - <s_i,e_u> + <s_i,c_u>: avg= -0.28125, std= 273454  [-4954579,4954578] (11672290)
-c_v: avg= 0.125, std= 2.42129e+06  [-4193792,4193792] (8387585)
PKE: <e,r> + y_i - <s_i,e_u> - c_v + <s_i,c_u>: avg= -0.15625, std= 2.43668e+06  [-9052056,9052056] (23344580)
Total sigma/q= 0.072628
tail per coefficient = 9.881679e-56  2^-182.723

=== mmKyberPKE Summary N= 1024
mmKyberPKE-128: 1.548895e-45  2^-148.856
mmKyberPKE-192: 3.801779e-65  2^-213.999
mmKyberPKE-256: 2.590423e-50  2^-164.723
```
