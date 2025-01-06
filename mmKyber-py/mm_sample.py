#   mm_sample.py
#   === Uniform and Discrete Gaussian samplers

import math
from mm_ring import MM_D, MM_Q

def poly_unif(xof, q=MM_Q, d=MM_D):
    """ Sample an uniform integer polynomial in [0,q-1]. """
    bits_q  = int(q - 1).bit_length()
    byte_q  = (bits_q + 7) // 8
    mask_q  = (1 << bits_q) - 1
    r   = []
    while len(r) < d:
        x   = int.from_bytes(xof.read(byte_q), byteorder='little') & mask_q
        if x < q:
            r += [x]
    return r

def sample_tri(xof, d=MM_D):
    """ Sample a polynomial in ternary uniform set {-1,0,1}."""
    b = b''
    i = 0
    while i < d:
        x = xof.read(1)[0]
        if x < 243:
            b += bytes([x])
            i += 5
    return b

def sample_bin(xof, d=MM_D):
    """ Sample a polynomial in the binary set {0,1}."""
    return xof.read((d + 7) // 8);

def poly_tri(s, d=MM_D):
    """ Decode a ternary polynomial from bytes in [0,243). """
    v = []
    for x in s:
        for i in range(5):
            v += [ (x % 3) - 1 ]
            x //= 3
    return v[:d]

def poly_bin(s, d=MM_D):
    """ Decode a binary polynomial from bytes. """
    v = []
    for x in s:
        v += [ (x >> i) & 1  for  i in range(8) ]
    return v[:d]


"""
=== IMPORANT NOTICE: This sampler is approximate and not constant-time;
    it is a placeholder implementation. A more appropriate Discrete Gaussian
    sampler is required in practice.

Firstly we note that "Gaussian width" is related to standard distribution
by sigma = gw / sqrt(2*Pi).

Rounded Gaussian sampler uses narrower sigma' = sqrt(sigma^2 - 1/12)
in the continuous sampler to get the desired distribution sigma.
For intuition, note that 1/12 is the variance of [-0.5, +0.5]
uniform distribution; models the rounding effect when sigma is larger.
Sect A.3 in https://eprint.iacr.org/2017/1025.pdf
Uses S. Janson 2006: http://dx.doi.org/10.1214/009117906000000232
For a Renyi argument, see Appendix E in https://eprint.iacr.org/2024/184.pdf

    orig    cs2 =   -2 * sigma'^2
    substitute      -2 * (sigma^2 - 1/12)
    simplify        1/6 - 2*sigma^2
    simplify        1/6 - 2*(gw/sqrt(2*Pi))^2
                    1/6 - gw^2/Pi

"""

def poly_gauss(xof, gw, d=MM_D):
    """ Discrete gaussian with width gw = sigma/sqrt(2*Pi). """
    cs2 =   1.0/6.0 - (gw*gw)/math.pi       #   adjust sigma
    d63 =   0.5**63                         #   scaling
    r   =   []
    while len(r) < d:
        x   = d63 * int.from_bytes(xof.read(8), byteorder='little') - 1.0
        y   = d63 * int.from_bytes(xof.read(8), byteorder='little') - 1.0
        w   = x*x + y*y
        if w > 0.0 and w <= 1.0:
            w   = math.sqrt( cs2 * math.log(w) / w )
            r   += [ round(x * w) ]
            r   += [ round(y * w) ]
    return r

