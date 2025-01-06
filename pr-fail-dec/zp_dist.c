//  zp_dist.c
//  === Discrete distributions

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

#include <mpfr.h>
#include <math.h>

//  discrete distributions
#include "crt_ntt.h"
#include "zp_dist.h"

//  how much space to actually allocate for product
#define DIST_PREC(x) (2 * (x) + 32)

//  initialize the system

int dist_init(size_t eps)
{
    (void) eps;
    crt_ntt_init();

    return 0;
}

//  Return a new distribution with "n" elements, precision "eps" bits.

dist_t *dist_alloc(size_t n, size_t eps)
{
    size_t i, prec;
    mpz_t *v;
    dist_t *d;

    d = calloc(1, sizeof(dist_t));
    assert(d != NULL);
    v = calloc(n, sizeof(mpz_t));
    assert(v != NULL);
    d->x = 0;
    d->n = 0;
    d->v = v;
    d->eps = eps;
    d->v_sz = n;
    prec = DIST_PREC(eps);
    for (i = 0; i < n; i++) {
        mpz_init2(v[i], prec);
    }

    return d;
}

//  Make sure distribution "d" has space for "sz" values, set to 0.

dist_t *dist_grow(dist_t *dr, size_t sz, size_t eps)
{
    size_t i, prec;

    //  newly allocated vector
    if (dr == NULL) {
        return dist_alloc(sz, eps);
    }

    //  already has sufficient space?
    if (dr->v_sz >= sz) {
        return dr;
    }

    //  at least double it (to reduce the number of reallocations)
    if (sz < dr->v_sz * 2) {
        sz = dr->v_sz * 2;
    }

    //  grow the vector
    dr->v = realloc(dr->v, sizeof(mpz_t) * sz);
    assert(dr->v != NULL);

    prec = DIST_PREC(dr->eps);          //  intialize new elements
    for (i = dr->v_sz; i < sz; i++) {
        mpz_init2(dr->v[i], prec);
    }
    dr->v_sz = sz;

    return dr;
}

//  zeroize a vector

dist_t *dist_zero(dist_t *dr, size_t n, size_t eps)
{
    size_t i;

    dr = dist_grow(dr, n, eps);

    for (i = 0; i < n; i++) {
        mpz_set_ui(dr->v[i], 0);
    }

    return dr;
}

//  set dr to da

dist_t *dist_copy(dist_t *dr, const dist_t *da)
{
    size_t i;

    dist_grow(dr, da->n, da->eps);
    dr->x = da->x;
    dr->n = da->n;
    for (i = 0; i < da->n; i++) {
        mpz_set(dr->v[i], da->v[i]);
    }
    return dr;
}

//  make a new copy

dist_t *dist_dup(const dist_t *da)
{
    return dist_copy(NULL, da);
}

//  Free distribution "d".

void dist_clear(dist_t *d)
{
    size_t i;

    if (d == NULL) {
        return;
    }
    if (d->v != NULL) {
        for (i = 0; i < d->v_sz; i++) {
            mpz_clear(d->v[i]);
        }
        free(d->v);
    }
    free(d);
}

//  print information about distribution "d", labeled "lab"
//  returns approximate variance

double dist_print(const dist_t *d, const char *lab)
{
    size_t i;
    int64_t x;
    double sum, avg, var, std;
    mpz_t az, bz, cz, xz;

    mpz_inits(az, bz, cz, xz, NULL);

    //  compute statistics
    for (i = 0; i < d->n; i++) {

        x = d->x + ((int64_t) i);           //  coord

        mpz_add(    az, az, d->v[i]);       //  mass:   a += y

        mpz_mul_si( xz, d->v[i], x);        //  avg:    b += x * y
        mpz_add(    bz, bz, xz);

        mpz_mul_si( xz, xz, x);             //  var:    c += x^2 * y
        mpz_add(    cz, cz, xz);
    }
    sum = mpz_get_d(az);
    avg = mpz_get_d(bz);
    var = mpz_get_d(cz);

    if (sum > 0.0) {
        avg = avg / sum;
        var = (var / sum) - (avg * avg);
        std = sqrt(var);
    } else {
        avg = 0.0;
        var = 0.0;
        std = 0.0;
    }

    mpz_clears(az, bz, cz, xz, NULL);

    printf("%s: avg= %g, std= %g  [%ld,%ld] (%zu)\n",
        lab, avg, std, d->x, d->x + ((int64_t) d->n) - 1, d->v_sz);
    fflush(stdout);

    return var;
}

//  dumps internal vector as a polynomial (suitable for gp/pari)

void vecz_dump(mpz_t *v, size_t n, const char *lab)
{
    size_t i;
    char buf[1000] = "0";

    if (n > 0) {
        mpz_get_str(buf, 10, v[0]);
    }
    printf("%s= %s", lab, buf);

    for (i = 1; i < n; i++) {
        if (mpz_sgn(v[i]) > 0) {
            mpz_get_str(buf, 10, v[i]);
            printf(" + %s*x^%zu", buf, i);
        }
    }
    printf(";\n");
}

//  Scale to sum to 2**eps, prune. If lab != NULL, report size.
size_t dist_sum1(dist_t *d, const char *lab)
{
    size_t i, lo, hi, n, drop;
    mpz_t az, xz;
    double sum;

    mpz_inits(az, xz, NULL);

    //  compute statistics
    for (i = 0; i < d->n; i++) {
        mpz_add(az, az, d->v[i]);           //  a += y
    }
    sum = mpz_get_d(az);

    //  bit drop
    drop = mpz_sizeinbase(az, 2);
    if (drop > d->eps) {
        drop = drop - d->eps;
    } else {
        drop = 0;
    }

    //  scaling factor
    mpz_set_si(xz, 0);
    mpz_setbit(xz, 2 * d->eps + drop);
    mpz_fdiv_q(az, xz, az);

    lo  = d->n;
    hi  = 0;
    for (i = 0; i < d->n; i++) {

        mpz_fdiv_q_2exp(xz, d->v[i], drop);
        mpz_mul(xz, xz, az);
        mpz_fdiv_q_2exp(d->v[i], xz, d->eps);

        if (mpz_sgn(d->v[i]) > 0) {
            if (lo > i)
                lo = i;
            hi = i;
        }
    }
    if (lo >= d->n) {
        lo = 0;         //  this actually means everything vanished!
    }
    n = hi - lo + 1;

    //  print statistics before adjustment
    if (lab != NULL) {
        printf("%s: sum= 2^%.15g, eps= %zu, trim [%ld,%ld]->[%ld,%ld] (%zu)\n",
            lab, log(sum)/log(2), d->eps, d->x, d->x + ((int64_t) d->n) - 1,
            d->x + ((int64_t) lo), d->x + ((int64_t) hi), d->v_sz);
    }

    //  pruning
    if (lo > 0) {
        for (i = 0; i < n; i++) {
            mpz_set(d->v[i], d->v[i + lo]);
        }
        for (i = n; i < d->n; i++) {
            mpz_set_ui(d->v[i], 0);
        }
    }
    d->n = n;
    d->x += lo;

    mpz_clears(az, xz, NULL);

    return n;
}


//  create a discrete gaussian distribution with given parameters
//  width=false: use std. devation. true: scale sigma by 1/sqrt(2*Pi)

dist_t *dist_gauss( dist_t *dr, const char *mu,
                    const char *sig, bool width,
                    size_t eps)
{
    int64_t x, lox, hix;
    size_t  i, n;
    mpfr_t muf, sigf, xf;

    mpfr_inits2(eps, xf, sigf, muf, NULL);

    //  parse (mu, sig) strings
    if (mpfr_set_str(muf, mu, 0, MPFR_RNDN) ||
        mpfr_set_str(sigf, sig, 0, MPFR_RNDN)) {
        fprintf(stderr, "dist_gauss() invalid inputs.\n");
        exit(1);
    }

    //  if width=true scale sigma by 1/sqrt(2*Pi)
    if (width) {
        mpfr_const_pi(xf, MPFR_RNDN);
        mpfr_mul_si(xf, xf, 2, MPFR_RNDN);
        mpfr_sqrt(xf, xf, MPFR_RNDN);
        mpfr_div(sigf, sigf, xf, MPFR_RNDN);
    }

    //  estimate tailcut
    //  solve: exp(-((x-mu)/sig)^2/2)/(sig*sqrt(2*Pi)) = y = 2^-eps
    //  for x: x = mu + sig*sqrt(-2*log(y*(sig*sqrt(2*Pi))))
    //           = mu + sig*sqrt(2)*sqrt(log(2)*eps - log(sig*sqrt(2*Pi)))
    double mud  = mpfr_get_d(muf, MPFR_RNDN);
    double sigd = mpfr_get_d(sigf, MPFR_RNDN);
    double tau  = sigd * M_SQRT2 * sqrt(M_LN2 * ((double) eps) -
                                        log(sigd*sqrt(2.0 * M_PI)));
    lox = floor(mud - tau);
    hix = ceil(mud + tau);
    n   = hix - lox + 1;

    //  allocate distribution
    dr  = dist_grow(dr, n, eps);
    dr->x = lox;
    dr->n = n;

    for (i = 0; i < n; i++) {

        x = lox + ((int64_t) i);

        //  compute exp(-((x-mu)/sig)^2/2)
        mpfr_set_si(xf, x, MPFR_RNDN);
        mpfr_sub(xf, xf, muf, MPFR_RNDN);
        mpfr_div(xf, xf, sigf, MPFR_RNDN);
        mpfr_mul(xf, xf, xf, MPFR_RNDN);
        mpfr_div_si(xf, xf, -2, MPFR_RNDN);
        mpfr_exp(xf, xf, MPFR_RNDN);

        //  scale by eps and store
        mpfr_mul_2ui(xf, xf, eps, MPFR_RNDN);
        mpfr_get_z(dr->v[i], xf, MPFR_RNDN);
    }

    mpfr_clears(xf, muf, sigf, NULL);

    //  normalize
    dist_sum1(  dr, NULL );

    return dr;
}

//  create a uniform distribution in range [a, b]

dist_t *dist_unif(dist_t *dr, int64_t a, int64_t b, size_t eps)
{
    size_t i, n;

    assert(a <= b);
    n = b - a + 1;
    dr = dist_grow(dr, n, eps);

    //  uniform
    dr->x = a;
    dr->n = n;

    for (i = 0; i < n; i++) {
        mpz_set_ui(dr->v[i], 1);
    }

    //  scale it
    dist_sum1(dr, NULL);

    return dr;
}

//  centered binomial distribution [-eta, eta]

dist_t *dist_cbd(dist_t *dr, int64_t eta, size_t eps)
{
    size_t n;
    int64_t x, w;

    n = 2 * eta + 1;
    dr = dist_zero(dr, n, eps);
    dr->n = n;
    dr->x = -eta;

    for (x = 0; x < (1 << (2 * eta)); x++) {
        w = __builtin_popcountl(x);
        mpz_add_ui(dr->v[w], dr->v[w], 1);
    }

    //  scale it
    dist_sum1(dr, NULL);

    return dr;
}

//  normalize value to [0,q-1]

static inline int64_t unsigned_q(int64_t x, int64_t q)
{
    x %= q;
    if (x < 0) {
        x += q;
    }

    return x;
}

//  normalize to [-ceil(q/2),floor(q/2)]

static inline int64_t signed_q(int64_t x, int64_t q)
{
    x %= q;
    if (x < 0) {
        x += q;
    }
    if (x >= (q/2)) {
        x -= q;
    }

    return x;
}

//  compress from (mod q) to (mod 2^d)

static inline int64_t compress_q(int64_t x, int64_t d, int64_t q)
{
    //  normalize to [0,q-1]
    x = unsigned_q(x, q);

    //  scale, add a rounding constant
    x = (x << d) + (q / 2l);

    //  divide by Q, normalize to [0, 2^d-1]
    x /= q;
    x &= (1l << d) - 1l;

    return x;
}

//  decompress from (mod 2^d) to (mod q)

static inline int64_t decompress_q(int64_t x, int64_t d, int64_t q)
{
    //  round up and add rounding constant
    x   &= (1l << d) - 1l;
    x   = (x * q) + (1l << (d - 1l));

    //  divide, normalize
    x   >>= d;

    //  normalize to [0,q-1]
    x = unsigned_q(x, q);

    return x;
}

//  Psi, the roudning distribution: r - Decompress_d(Compress_d(r))

dist_t *dist_round(dist_t *dr, int64_t d, int64_t q, size_t eps)
{
    size_t n;
    int64_t r, x, b;

    b   = (q >> (d + 1)) + 1;   //  allocate [-b,+b]
    n   = 2 * b + 1;

    dr  = dist_zero(dr, n, eps);
    dr->x = -b;
    dr->n = n;

    //  we could of course make this quicker with congruence classes..
    for (r = 0; r < q; r++) {
        x = r - decompress_q(compress_q(r, d, q), d, q);
        x = signed_q(x, q);
        assert(x >= -b && x <= b);
        x += b;
        mpz_add_ui(dr->v[x], dr->v[x], 1);
    }

    //  normalize
    dist_sum1(dr, NULL);

    return dr;
}

//  Double-and-round distribution: c_i - [dbl(c_i]_{2^{d_u}}

dist_t *dist_dbl_r(dist_t *dr, int64_t d, int64_t q, size_t eps)
{
    size_t n;
    int64_t r, x, b;

    b   = (q >> d) + 1; //  allocate [-b,+b]
    n   = 2 * b + 1;

    dr  = dist_zero(dr, n, eps);
    dr->x = -b;
    dr->n = n;

    for (r = 0; r < q; r++) {

        //  model the dbl() randomness here
        x = 2 * r;
        x = x - decompress_q(compress_q(x, d, 2 * q), d, 2 * q);
        x = signed_q(x, 2 * q);
        assert(x >= -b && x <= b);
        x += b;
        mpz_add_ui(dr->v[x], dr->v[x], 2);      //  2*x:    Pr. 2/4

        x = 2 * r - 1;
        x = x - decompress_q(compress_q(x, d, 2 * q), d, 2 * q);
        x = signed_q(x, 2 * q);
        assert(x >= -b && x <= b);
        x += b;
        mpz_add_ui(dr->v[x], dr->v[x], 1);      //  2*x-1:  Pr. 1/4

        x = 2 * r + 1;
        x = x - decompress_q(compress_q(x, d, 2 * q), d, 2 * q);
        x = signed_q(x, 2 * q);

        assert(x >= -b && x <= b);
        x += b;
        mpz_add_ui(dr->v[x], dr->v[x], 1);      //  2*x+1:  Pr. 1/4
    }

    //  normalize
    dist_sum1(dr, NULL);

    return dr;
}

//  negate a distribution: r = -a

dist_t *dist_neg(dist_t *dr, const dist_t *da)
{
    size_t i, n, eps;
    dist_t *rr = NULL;

    //  allocate result
    n   = da->n;
    eps = da->eps;
    if (dr == NULL || dr == da) {
        rr = dist_alloc(n, eps);
    } else {
        rr = dist_grow(dr, n, eps);         //  not zeroed!
    }

    for (i = 0; i < n; i++) {
        mpz_set(rr->v[i], da->v[n - 1 - i]);
    }
    rr->n = n;
    rr->x = -(da->x + ((int64_t) n) - 1);   // new first = -last

    //  free variables
    if (dr == NULL) {
        return rr;
    }
    if (dr != rr) {
        dist_copy(dr, rr);
        dist_clear(rr);
    }
    return dr;
}

//  Add (convolute) two distributions, write to dr

dist_t *dist_add(dist_t *dr, const dist_t *da, const dist_t *db)
{
    size_t  n, eps;
    dist_t *rr = NULL;

    //  length is sum; precision (if needed); pick larger
    n   = da->n + db->n - 1;
    eps = da->eps > db->eps ? da->eps : db->eps;

    if (da->n == 0 || db->n == 0) {
        dr = dist_zero(dr, 1, eps);
        dr->n = 0;
        return dr;
    }

    //  allocate result
    if (dr == NULL || dr == da || dr == db) {
        rr = dist_alloc(n, eps);
    } else {
        rr = dist_zero(dr, n, eps);
    }

    rr->x = da->x + db->x;
    rr->n = n;

    //  NTT convolution
    crt_ntt_convol(rr->v,
                da->v, da->n,
                db->v, db->n,
                DIST_PREC(dr->eps));

    //  normalize
    dist_sum1(rr, NULL);

    //  free variables
    if (dr == NULL) {
        return rr;
    }
    if (dr != rr) {
        dist_copy(dr, rr);
        dist_clear(rr);
    }
    return dr;
}

//  Signed multiply of two distributions

dist_t *dist_mul(dist_t *dr, const dist_t *da, const dist_t *db)
{
    int64_t x, y, t;
    int64_t ax, ay, bx, by;
    size_t  i, j, n, eps;
    dist_t  *rr = NULL;
    mpz_t   rz;

    mpz_inits(rz, NULL);

    //  get bounds of the product.
    ax  = da->x;    //  da in [ax,ay]
    ay  = ax + ((int64_t) da->n) - 1;
    bx  = db->x;    //  db in [bx,by]
    by  = bx + ((int64_t) db->n) - 1;


    //  d in [x,y], start with a possible upper or lower bound
    //  we may have x>0 or y<0 situations, hence up/low bound can be any
    t   = ax * bx;
    x   = t;
    y   = t;

    t   = ax * by;
    x   = t < x ? t : x;
    y   = t > y ? t : y;

    t   = ay * bx;
    x   = t < x ? t : x;
    y   = t > y ? t : y;

    t   = ay * by;
    x   = t < x ? t : x;
    y   = t > y ? t : y;

    //  range
    n   = y - x + 1;

    //  big the bigger precision
    eps = da->eps > db->eps ? da->eps : db->eps;

    //  allocate result
    if (dr == NULL || dr == da || dr == db) {
        rr = dist_alloc(n, eps);
    } else {
        rr = dist_zero(dr, n, eps);
    }
    rr->x = x;
    rr->n = n;

    //  multiply probability masses
    for (i = 0; i < da->n; i++) {
        if (mpz_sgn(da->v[i]) == 0)
            continue;
        ax = da->x + ((int64_t) i);

        for (j = 0; j < db->n; j++) {
            if (mpz_sgn(db->v[j]) == 0)
                continue;
            bx = db->x + ((int64_t) j);

            t = ax * bx;
            assert(t >= x && t <= y);

            t -= x;
            mpz_mul(rz, da->v[i], db->v[j]);
            mpz_add(rr->v[t], rr->v[t], rz);
        }
    }

    //  free variables
    mpz_clears(rz, NULL);

    dist_sum1(rr, NULL);

    if (dr == NULL) {
        return rr;
    }
    if (dr != rr) {
        dist_copy(dr, rr);
        dist_clear(rr);
    }

    return dr;
}

//  Add a distribution to itself x-1 times ("scalar multiplication" by x)

dist_t *dist_scmul(dist_t *dr, const dist_t *da, uint64_t x)
{
    size_t n, eps;
    dist_t *dt = NULL, *rr = NULL;
    int i;

    //  if x == 0 return r := 0 distribution
    if (x == 0) {
        return dist_unif(dr, 0, 0, da->eps);
    }

    //  need to allocate result?
    n   = da->n;                    //  (will be grown later)
    eps = da->eps;
    if (dr == NULL || dr == da) {
        rr = dist_alloc(n, eps);
    } else {
        rr = dist_zero(dr, n, eps);
    }

    //  find the top bit
    i = 63;
    while (i > 0 && (x >> i) == 0) {
        i--;
    }

    dist_copy(rr, da);              //  r := a   set top bit
    dt = dist_alloc(n, eps);        //  temporary
    i--;
    while (i >= 0) {

        dt = dist_add(dt, rr, rr);  //  double t := 2*r

        if ((x >> i) & 1) {
            dist_add(rr, dt, da);   //  if 1, r := 2*r + a
        } else {
            dist_copy(rr, dt);      //  if 0, r := 2*r
        }
        i--;
    }

    //  free variables
    dist_clear(dt);
    if (dr == NULL) {
        return rr;
    }
    if (dr != rr) {
        dist_copy(dr, rr);
        dist_clear(rr);
    }
    return dr;
}

//  multiply a distribution by a constant, spreading it by factor c
dist_t *dist_spread(dist_t *dr, const dist_t *da, int64_t c)
{
    size_t i, n, eps;
    dist_t *rr = NULL;

    //  note: we could handle negative c with dist_neg() here
    assert(c > 0);

    //  trivial distributions
    if (da->n <= 1) {
        rr = dist_copy(dr, da);
        rr->x *= c;
        return rr;
    }

    //  allocate result
    n   = c * (da->n - 1) + 1;
    eps = da->eps;
    if (dr == NULL || dr == da) {
        rr = dist_alloc(n, eps);
    } else {
        rr = dist_zero(dr, n, eps);
    }

    //  spread it
    for (i = 0; i < da->n; i++) {
        mpz_set(rr->v[c * i], da->v[i]);
    }
    rr->n = n;
    rr->x = c * da->x;

    //  free variables
    if (dr == NULL) {
        return rr;
    }
    if (dr != rr) {
        dist_copy(dr, rr);
        dist_clear(rr);
    }
    return dr;
}
//  return tail -- probability mass outside (lo, hi)

double dist_tail(dist_t *d, int64_t lo, int64_t hi)
{
    size_t i;
    int64_t x;
    mpz_t az, bz;
    double mass, tail;

    mpz_inits(az, bz, NULL);

    //  compute statistics
    for (i = 0; i < d->n; i++) {

        x = d->x + ((int64_t) i);           //  coord

        mpz_add(az, az, d->v[i]);           //  total: a += y
        if (x <= lo || x >= hi) {
            mpz_add(bz, bz, d->v[i]);       //  tail: b += y
        }
    }

    mass = mpz_get_d(az);                   //  can lower precision now
    tail = mpz_get_d(bz);
    mpz_clears(az, bz, NULL);

    return tail / mass;
}

