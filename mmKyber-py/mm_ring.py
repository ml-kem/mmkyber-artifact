#   mm_ring.py
#   === General Polynomial Ring Arithmetic (+vector built on top.).

MM_D    = 256           #   base polynomial x^d + 1
MM_Q    = 33550337      #   2**25 - 2**12 + 1
MM_H    = 8433925       #   2d'th root of unity
MM_DI   = 33419281      #   d^-1  (mod q)

#   These constants can be calculated with _calc_w() in test_ntt.py

MM_W = [
    1,          12759331,   13682589,   27193631,   17371284,   22604314,
    31280465,   11836922,   30292897,   6209026,    878186,     24953317,
    15930566,   23460663,   14378272,   30651277,   255445,     26269093,
    9039793,    13996293,   6519423,    11790682,   23021331,   25518839,
    18697474,   6015232,    10669588,   5084772,    29506803,   13663747,
    1648639,    6206901,    27306318,   30908010,   10619470,   2356934,
    1754428,    5915876,    18882277,   17427328,   28049154,   16906878,
    23459409,   20929446,   20460260,   31448620,   15531894,   3346823,
    13137862,   30009588,   11566352,   6208165,    28009151,   6665666,
    19049460,   30235441,   6173810,    10320385,   28839087,   9031846,
    33168177,   32943946,   16010558,   33064138,   14921815,   18439066,
    29281037,   22476021,   8236665,    18089846,   32009659,   27740981,
    5983912,    10669950,   4851793,    27212911,   30863996,   16140428,
    5433838,    30728508,   15695768,   1852603,    15916022,   18664546,
    6156981,    5695786,    20511637,   27562764,   7047120,    23100544,
    16814105,   7076354,    26216253,   24266867,   2205546,    6881540,
    7702823,    7828558,    2301665,    21999231,   22736731,   7091097,
    375524,     10736463,   8081092,    16376788,   25272464,   29830107,
    6160655,    13657091,   31938797,   24059472,   21007196,   31711762,
    12710337,   12766306,   18311551,   2578535,    5314697,    3492970,
    22961341,   10640467,   2271277,    32693512,   29959490,   4468561,
    2799690,    20442369,   8433925,    10991503,   2004823,    14210459,
    25093382,   9017350,    26285948,   608390,     24303820,   20006666,
    16014267,   17897310,   9536185,    11211185,   23765701,   24484990,
    2631507,    30981653,   9667267,    11987540,   19329455,   5033678,
    17291365,   5022566,    740072,     14162508,   5393742,    8131308,
    15009103,   17936342,   30213143,   13795999,   27520735,   19702709,
    14636455,   32517494,   9042790,    6327120,    926175,     5288089,
    28452340,   26712223,   3425390,    29857897,   11269975,   26648996,
    21012388,   16242739,   30738443,   5346461,    26792869,   5270033,
    28339437,   10784099,   23346688,   12226311,   32037327,   2714838,
    5959590,    18838618,   4996916,    22956257,   25888389,   23387739,
    7265307,    3389855,    3773314,    27395912,   32987471,   25085911,
    15718276,   31540368,   5434709,    5393925,    23305812,   1801020,
    9856327,    21876100,   29624282,   5153810,    15905123,   20862842,
    6563057,    18147358,   15439012,   23269069,   18432345,   17826243,
    23396119,   5931709,    13597375,   19332956,   31511024,   29784117,
    29104466,   33321907,   4172010,    18230326,   12748610,   25486634,
    13996356,   23140309,   22989237,   9412810,    28552842,   4015571,
    3322449,    158302,     1168074,    11445343,   3077907,    22728237,
    26189982,   25299133,   5220545,    23877617,   7496015,    10108160,
    2408670,    3248671,    5212575,    23080994,   13660053,   9298305,
    15515989,   12175781,   17356357,   29333626    ]

### Internal functions: Polynomial arithmetic

def ntt(f, d=MM_D, w=MM_W, q=MM_Q):
    """ Forward NTT (negacyclic - x^d+1.) Note: Transforms f in place."""

    l = d // 2
    wi = 0
    while l > 0:
        for i in range(0, d, 2 * l):
            wi += 1
            z = w[wi]
            for j in range(i, i + l):
                x = f[j]
                y = (f[j + l] * z) % q
                f[j] = (x + y) % q
                f[j + l] = (x - y) % q
        l >>= 1
    return f


def intt(f, d=MM_D, w=MM_W, di=MM_DI, q=MM_Q):
    """ Inverse NTT (negacyclic - x^n+1.) Note: Transforms f in place."""

    wi = d
    l = 1
    while l < d:
        for i in range(0, d, 2 * l):
            wi -= 1
            z = w[wi]
            for j in range(i, i + l):
                x = f[j]
                y = f[j + l]
                f[j] = (x + y) % q
                f[j + l] = (z * (y - x)) % q
        l <<= 1
    #   normalize: ni = n^-1  (mod q)
    for i in range(d):
        f[i] = (di * f[i]) % q
    return f

def vec_ntt(v):
    """ NTT on a vector-like object of polynomials (copying)."""
    return [ ntt(vi.copy()) for vi in v ]

def vec_intt(v):
    """ NTT^-1 on a vector-like object of polynomials (copying)."""
    return [ intt(vi.copy()) for vi in v ]


def mul_ntt(f, g, q=MM_Q):
        """ Multiplication of two polynomials (NTT domain.)"""
        return [ (fi * gi) % q for fi,gi in zip(f,g) ]

"""

#   Kyber-style middle-layer method (if we can only do 7 layers of NTT.)
def mul_ntt(f, g, d=MM_D, w=MM_W, q=MM_Q):
    fg  = []
    for i in range(0, d, 4):
        z2  = w[64 + i//4]

        f0  = f[i]
        f1  = f[i + 1]
        g0  = g[i]
        g1  = g[i + 1]
        fg0 = (g1*f1*z2 + g0*f0) % q
        fg1 = (g1*f0 + g0*f1) % q
        fg  += [ fg0, fg1 ]

        f0  = f[i + 2]
        f1  = f[i + 3]
        g0  = g[i + 2]
        g1  = g[i + 3]
        fg0 = (-g1*f1*z2 + g0*f0) % q
        fg1 = (g1*f0 + g0*f1) % q
        fg  += [ fg0, fg1 ]
    return fg
"""

def poly_zero(d=MM_D):
    """ Zero polynomial. """
    return [ 0 for _ in range(d) ]

def poly_add(f, g, q=MM_Q):
    """ Add polynomials: return f + g (mod q)."""
    return [ (fi + gi) % q for fi,gi in zip(f,g) ]

def vec_add(t, u, q=MM_Q):
    """ Add vectors of polynomials: return t + u."""
    return [ poly_add(ti, ui, q=q) for ti,ui in zip(t, u) ]

def poly_sub(f, g, q=MM_Q):
    """ Subtract polynomials: return f - g (mod q)."""
    return [ (fi - gi) % q for fi,gi in zip(f,g) ]

def vec_sub(t, u, q=MM_Q):
    """ Subtract vectors of polynomials: return t - u."""
    return [ poly_sub(ti, ui, q=q) for ti,ui in zip(t, u) ]

def poly_scale(c, f, q=MM_Q):
    """ Multiplication of a polynomials with scalar: c * f  (mod q).)"""
    return [ (c * fi) % q for fi in f ]

def vec_scale(c, t, q=MM_Q):
    """ Multiplication of a vector of polys with scalar: c * t (mod q).)"""
    return [ poly_scale(c, ti, q=q) for ti in t ]

def vec_mul_ntt(f, v, q=MM_Q):
    """ Multiplication of polynomial f by vector v (NTT domain.)"""
    return [ mul_ntt(f, vi, q=q) for vi in v ]

def poly_signed(f, q=MM_Q):
    """ Normalize a polynomial coefficients to [-q/2, q/2]. """
    return [ (x + q//2) % q - q//2  for x in f ]

def vec_signed(t, q=MM_Q):
    """ Normalize vector coefficients to [-q/2, q/2]. """
    return [ poly_signed(ti, q=q) for ti in t ]

def poly_mul1_add(r, f, g):
    """ Compute r += f*g  mod x^d+1, no reduction. """
    d = len(r)
    for i in range(d):
        x = g[i]
        for j in range(i, d):
            r[j] += x * f[j - i]
        for j in range(i):
            r[j] -= x * f[d + j - i]

