#   test_misc.py
#   === some unit tests for mmKyber

from Crypto.Hash import SHAKE128
from random import randrange
from mm_ring import *
from mm_serial import *
from mm_sample import *

if (__name__ == "__main__"):

    def _modexp(x, e, n):
        """(TESTING) Modular exponentiation: Compute x**e (mod n)."""
        y = 1
        while e > 0:
            if e & 1 == 1:
                y = (y * x) % n
            x = (x * x) % n
            e >>= 1
        return y

    def _bitrev(x, l):
        """(TESTING) Return x with bits 0,1,..(l-1) in reverse order."""
        y = 0
        for i in range(l):
            y |= ((x >> i) & 1) << (l - i - 1)
        return y

    #   x=3 is the smallest value that gives h with right order:
    #   h   = Mod(x,q)^((q-1)/(2*d)) = 8433925
    #   check  h^d == -1

    def _calc_w():
        """(TESTING) Re-generate the NTT "tweak" table."""
        q   = MM_Q
        lgd = 8                 #   log2(n)
        d   = 2**lgd            #   length of the transform
        h   = MM_H          #   Generates a subgroup of order 2*d
        w   = []
        for i in range(d):
            j = _bitrev(i, lgd)
            x = (_modexp(h, j, q)) % q
            w.append(x)
        return w

    #   prettyprint the table
    def _print_w(r=1):
        w = _calc_w()
        s = '\t'
        for i in range(0, len(w)):
            x = (r * w[i]) % MM_Q
            n = str(x) + ','
            s += f'{n:12}'
            if i % 6 == 5:
                s += '\n\t'
        print(s)

    def _rand_poly(n=MM_D, q=MM_Q):
        """(TESTING) Random polynomial."""
        return [ randrange(q) for _ in range(n) ]

    def _conv_slow(f, g, d=MM_D,q=MM_Q):
        """(TESTING) Slow negacyclic convolution h = f*g (mod x^d+1)."""
        h = [0] * d
        for i in range(d):
            for j in range(d):
                x = (f[i] * g[j]) % q
                k = i + j
                if k >= d:
                    k -= d
                    x = -x                  # x^n == -1 (mod x^n + 1)
                h[k] = (h[k] + x) % q

        return h

    def _conv_fast(f, g, d=MM_D, q=MM_Q):
        """(TESTING) Fast NTT negacyclic convolution h = f*g (mod x^d+1)."""
        ft = ntt(f.copy())
        gt = ntt(g.copy())
        ht = mul_ntt(ft, gt)
        return intt(ht)

    def _test_rec(q2 = 2*MM_Q):
        """ Ding-Peikert Reconciliation test."""
        fail = 0
        q8  = q2 // 8   # test set for E
        eset = [ -q8, -q8 + 1, -1, 0, 1, q8 - 2, q8 - 1 ]

        for v in range(0, q2, 1000):        #   step 1000 for quick check
            round2 = poly_compress([v], 1, q2)[0]
            cross2 = poly_cross([v], q2)[0]
            for e in eset:
                w   = (v + e) % q2
                rec = poly_rec([w], [cross2], q2)[0]
                if rec != round2:
                    print(f'e={e}  w= {w}  b={cross2}')
                    fail += 1
        return fail

    def _stat(sampler, xof, arg, n=0x10000):
        """ Basic statistics from a sampler. """
        s   = 0.0
        r   = 0.0
        lo  = +1E100
        hi  = -1E100
        i   = 0
        while i < n:
            v = sampler(xof, arg)
            for x in v:
                lo  = min(lo, x)
                hi  = max(hi, x)
                s   += x
                r   += x * x
            i += len(v)
        avg = s / i
        std = (r / i - avg * avg)**0.5
        return f'{arg} n={i} [{lo},{hi}] avg={avg} std={std}'

    def test_sample():
        xof = SHAKE128.new(b'0')
        print('poly_gauss():', _stat( poly_gauss, xof, 15.90, 0x100000))
        print('poly_gauss():', _stat( poly_gauss, xof, 368459.34, 0x100000))

    def test_compress(dx, test=0):
        n       = 256
        q       = MM_Q
        fail    = 0

        xof = SHAKE128.new(test.to_bytes(8))
        v0  = poly_unif(xof, MM_Q)
        sl  = (n * dx + 7) // 8
        c0  = poly_compress(v0, dx)
        s   = poly_serial(c0, dx)
        if len(s) != sl:
            print(f'poly_serial(): dx={dx} n={n} |{len(s)}| bytes ({sl})')
            fail += 1
        c1, si  = poly_deserial(s, dx, n)
        if len(s) != sl:
            print(f'poly_deserial() dx={dx} n={n}  index={si} ({sl})')
            fail += 1
        if c0 != c1:
            print(f'poly_deserial() != poly_serial()')
            fail += 1
        v1  = poly_decompress(c1, dx)
        #   difference
        vd  = [ (v0[i] - v1[i] + q//2) % q - q//2 for i in range(n) ]
        h   = ((q >> dx) + 1) >> 1
        if min(vd) < -h:
            print(f'poly_compress() dx={dx}  bad rounding {min(vd)} < {-h}')
            fail += 1
        if max(vd) > h:
            print(vd)
            print(f'poly_compress() dx={dx}  bad rounding {max(vd)} > {h}')
            fail += 1

        return fail

    #   ---------------
    #_print_w(1)
    #_print_w(2**32)
    #exit(0)

    print('test_compress():', end='')
    for i in range(1,26):
        print(f' {i}: {test_compress(i)==0}', end='')
    print()

    print('_test_rec():', _test_rec() == 0)

    print("_calc_w():", _calc_w() == MM_W)

    #   test convolutions
    for _ in range(1):
        f = _rand_poly()
        g = _rand_poly()
        print("_conv_fast():", _conv_slow(f, g) == _conv_fast(f, g))

    test_sample()
