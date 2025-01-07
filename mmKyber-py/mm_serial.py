#   mm_serial.py
#   === Serialization and the Ding-Peikert Reconciliation Mechanism.

from mm_ring import MM_D,MM_Q

def poly_dbl(v, xof, d=MM_D, q=MM_Q):
    """
    Randomly double the modulus.
    return 2v + e \\mod 2q, e \\in {-1, 0, 1} in probability {1/4, 1/2, 1/4}
    """
    q2 = 2*q
    r = []
    for i in range(0, d, 32):
        #   64-bit chunk from XOF
        bits = int.from_bytes(xof.read(8), byteorder='little')

        #   e = bit0 - bit1
        for j in range(32):
            bit0 = bits & 1
            bit1 = (bits >> 1) & 1
            r += [ ( 2 * v[i + j] + bit0 - bit1) % q2 ]
            bits >>= 2
    return r

def poly_compress(v, dx, q=MM_Q):
    """ Compress the polynomial to dx bits per coefficient. """
    t = 1 << dx
    h = q // 2
    return [ ((t * x + h) // q) % t  for x in v]

def poly_decompress(v, dx, q=MM_Q):
    """ Decompress the polynomial by scaling each coefficient. """
    h = 1 << (dx - 1)
    return [ ((q * x + h) >> dx) % q  for x in v ]

def poly_serial(p, dx):
    """ Serialize polynomial p to bytes, dx bits per coefficient. """
    m = (1 << dx) - 1
    s = b''
    l = 0
    t = 0
    for x in p:
        t |= (x & m) << l
        l += dx
        while l >= 8:
            s += bytes([t & 0xFF])
            t >>= 8
            l -= 8
    if l > 0:
        s += bytes([t & 0xFF])
    return s

def vec_serial(v, dx):
    """ Serialize a vector of polynomials, dx bits per coefficient."""
    s = b''
    for p in v:
        s += poly_serial(p, dx)
    return s

def poly_deserial(s, dx, d=256):
    """ Deserialize bytes as d-1 degree poly, dx bits per coefficient. """
    m = (1 << dx) - 1
    p = []
    l = 0
    t = 0
    i = 0
    while len(p) < d:
        while l < dx:
            t |= s[i] << l
            i += 1
            l += 8
        p += [ t & m ]
        l -= dx
        t >>= dx
    return p, i

def vec_deserial(s, dx, n, d=256):
    """ Deserialize n polynomials, degree d-1, dx bits per coefficient. """
    v = []
    i = 0
    for _ in range(n):
        p, l = poly_deserial(s[i:], dx, d)
        v += [ p ]
        i += l
    return v, i

def poly_cross(v, q2=2*MM_Q):
    """ Cross rounding. """
    return [ ( (4 * x) // q2 ) % 2 for x in v ]

def poly_rec(w, b, q2=2*MM_Q):
    """ Reconciliation function. Note: Code uses a simplified method. """

    #   Bounds for I0 + E
    i0l = ((q2 + 2) // 4) + (q2 // 8) - 1
    i0h = q2 - (q2 // 8) - 1

    #   Bounds for I1 + E
    i1l = (q2 // 8) - 1
    i1h = q2 - (q2 // 4) - (q2 // 8) - 1

    r = []
    for i in range(len(b)):
        x = w[i]
        if b[i] == 0:
            if x < i0l or x >= i0h:
                y = 0
            else:
                y = 1
        else:
            if x < i1l or x >= i1h:
                y = 0
            else:
                y = 1
        r += [y]

    return r

