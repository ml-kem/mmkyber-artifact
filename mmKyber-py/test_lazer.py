#   test_lazer.py
#   === Create a header file with the raw secret key parametrs (for testing).
#       note: Our running LaZer deployment doesn't need these header files.

from mm_ring import *
from mm_sample import *
from mm_serial import *

def lazer_h(iut, fn, pp, pk, sk, qi=MM_Q):
    """ Output secret key data in a format suitable for LaZer.  """
    a   = pp

    #   expand secret key s
    s   = [ iut.poly_nu(sk[ i*iut.sz_nu : (i+1)*iut.sz_nu ])
                for i in range(iut.m) ]

    p1  =   poly_zero() #   identity poly polynomial in ntt domain
    p1[0] = 1;
    ntt(p1);

    #   deserialize  public key
    t,_ = vec_deserial(pk, iut.logq, iut.n)

    #   convert public to normal domain
    t   = [ intt( mul_ntt(ti, p1) ) for ti in t ]

    #   b := A^T * s

    st  = vec_ntt(s)
    b   = [ poly_zero() for _ in range(iut.n)]
    for i in range(iut.n):
        for j in range(iut.m):
            b[i] = poly_add( b[i], mul_ntt(a[j][i], st[j]) )

    b   = vec_intt(b)

    #   e = t - A^T *s
    e   = vec_signed( vec_sub(t, b) )


    ds  = f'//  sec={iut.sec}, m={iut.m}, n={iut.m}, q_i={qi}\n\n'
    ds  += '#include <stdint.h>\n\n'

    #   convert from NTT, flatten A, reduce mod q_i
    ds  += 'static int64_t mm_a[] = { '
    for i in range(iut.n):
        for j in range(iut.m):
            ap = intt( mul_ntt(a[j][i], p1) )
            for x in ap:
                ds += f'{x % qi}, '
    ds += '};\n\n'

    #   flatten t, reduce mod q_i
    ds  += 'static int64_t mm_t[] = { '
    for ti in t:
        for x in ti:
            ds += f'{x % qi}, '
    ds += '};\n\n'


    #   flatten s and e, compute l2 norm

    ds  += 'static int64_t mm_se[] = { '
    l2  = 0
    for ei in (s + e):
        for x in ei:
            ds += f'{x}, '
            l2 += x**2
    ds += '};\n\n'

    #   write to file

    with open(fn, 'wb') as f:
        f.write(ds.encode('ascii'))

    return l2**0.5

