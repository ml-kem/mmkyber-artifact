#   mmKyber.py
#   === multi-message, multi-recipient public key encryption

from Crypto.Hash import SHAKE128,SHAKE256
from mm_ring import *
from mm_serial import *
from mm_sample import *

class mmKyber:

    #   initialize
    def __init__(self,  sec=128, N=2**10, q=MM_Q, d=256, m=4, n=4,
                        nu_bar=3, du=14, dv=2,
                        sigma0=15.90, sigma1=368459.34):

        self.sec    = sec                   #   security level {128, 192, 256}
        self.N      = N                     #   max number of recipients
        self.q      = q                     #   modulus
        self.logq   = (q-1).bit_length()    #   bit size of prime field
        self.d      = d                     #   ring polynomial x^d+1 degree
        self.m      = m                     #   matrix A rows
        self.n      = n                     #   matrix A columns
        self.nu_bar = nu_bar                #   range of secret uniform
        self.du     = du                    #   bit size of u
        self.dv     = dv                    #   bit size of v
        self.sigma0 = sigma0                #   Gaussian (r, e_u)
        self.sigma1 = sigma1                #   Gaussian y_i
        if self.nu_bar == 2:                #   binary secret
            self.sz_nu      = self.d // 8
            self.sample_nu  = sample_bin
            self.poly_nu    = poly_bin
        elif self.nu_bar == 3:              #   ternary secret
            self.sz_nu      = (self.d + 4) // 5
            self.sample_nu  = sample_tri
            self.poly_nu    = poly_tri
        else:
            print(f'Please provide a sampler for nu_bar={self.nu_bar}')

        #   A matrix always with SHAKE128, others depending on security level.
        self.xof_a  = SHAKE128
        if sec <= 128:
            self.xof = SHAKE128
        else:
            self.xof = SHAKE256

    def init_xof(self, seed):
        """ (Mainly a hook to dump random seed values.) """
        #print(seed.hex(), chr(seed[-1]))
        return self.xof.new(seed)


    def mmSetup(self, seed_a):
        """ mmSetup(1^lambda, N): Generate public parameter A from seed. """
        #   A <- U(R_q^m*n)  -- direct in NTT domain
        a   = []
        for i in range(self.m):
            v = []
            for j in range(self.n):
                xof = self.xof_a.new(seed_a + bytes([i,j]) + b'A')
                v += [ poly_unif(xof, self.q) ]
            a += [v]

        #   return pp := A
        pp  = a
        return pp


    def mmKGen(self, pp, seed_sk):
        """ mmKGen(pp): Generate an individual public key. """
        a   = pp

        #   (s, e) <- U(Snu^m) x U(Snu^n)
        xof = self.xof.new(seed_sk + b'S')
        sk  = [ self.sample_nu(xof)   for _ in range(self.m) ]
        s   = [ self.poly_nu(sk[i]) for i in range(self.m) ]

        xof = self.xof.new(seed_sk + b'E')
        e   = [ self.poly_nu(self.sample_nu(xof)) for _ in range(self.n) ]

        #   b := A^T * s + e
        s   = vec_ntt(s)
        b   = vec_ntt(e)
        for i in range(self.n):
            for j in range(self.m):
                b[i] = poly_add( b[i], mul_ntt(a[j][i], s[j]) )

        #   return ( pk := t, sk := (s, e) )
        pk  = vec_serial(b, self.logq)

        #   sk doesn't actully contain e, but "base-243" polynomials of s
        sk  = b''.join(sk)
        return pk, sk

    #   """ mmDecap(pp, sk, ct): Decapsulate individual secret. """
    """
    def mmDecap(self, pp, sk, ct):
        s   = [ self.poly_nu(sk[ i*self.sz_nu : (i+1)*self.sz_nu ])
                    for i in range(self.m) ]
        (u, ui) = ct
        u,_ = vec_deserial(u, self.du, self.m)
        ui,_ = poly_deserial(ui, 1)

        #   c' := [u mod 2^du](q)
        c   = [ poly_decompress(u[i], self.du) for i in range(self.m) ]
        c   = vec_ntt(c)

        #   w := [ <c',s> mod q ]_2q
        s   = vec_ntt(s)
        w   = poly_zero()
        for i in range(self.m):
            w = poly_add( w, mul_ntt(c[i], s[i]) )
        w   = intt(w)
        w2  = poly_scale(2, w, 2 * self.q)

        #   return K := mu<-rec(w, u)
        mu  = poly_rec(w2, ui)
        k   = poly_serial(mu, 1)
        return  k
    """

    def mmDecap(self, pp, sk, ct):
        """ mmDecap(pp, sk, ct): Decapsulate *without* NTT. """

        def _rec(w, b, q):
            """ Simplified reconciliation function."""
            t = (8 * w) // q
            t = (t + 2 * b + 1) % 8
            return t // 4

        #   deserialize
        s   = [ self.poly_nu(sk[ i*self.sz_nu : (i+1)*self.sz_nu ])
                    for i in range(self.m) ]
        (u, ui) = ct
        u,_ = vec_deserial(u, self.du, self.m)
        ui,_ = poly_deserial(ui, 1)

        #   w := <c',s> mod 2^u_i
        w   = poly_zero()
        for i in range(self.m):
            poly_mul1_add(w, u[i], s[i])

        #   return K := mu<-rec(w, u)
        mu  = []
        for i in range(self.d):
            mu += [ _rec(w[i], ui[i], 2**self.du) ]
        k   = poly_serial(mu, 1)
        return  k


    def mmEncap(self, pp, pkl, seed_e):
        """ mmEncap(pp, (pk_i) for i in [N]): Encapsulate to N recipients."""

        #   r := (r, e_u) <= D^n_sigma0 x D^M_sigma0
        xof = self.init_xof(seed_e + b'R')
        r_u = [ poly_gauss(xof, self.sigma0, self.d) for _ in range(self.n) ]
        xof = self.init_xof(seed_e + b'e')
        e_u = [ poly_gauss(xof, self.sigma0, self.d) for _ in range(self.m) ]
        r   = (r_u, e_u)

        #   ^ct <- mmEnc^i(pp; r)
        ctu = self.mmEnc_i(pp, r)

        self.mmEnc_i(pp, r)
        ctl = []
        kl  = []

        #   for i in [N]
        for i in range(len(pkl)):
            #   r_i := y_i <- D_sigma1
            xof = self.init_xof(seed_e +
                int(i).to_bytes(8, byteorder='little') + b'r')
            r_i = poly_gauss(xof, self.sigma1, self.d)

            #   (~ct_i, K_i) <- mmEncap^d(pp. pk_i; r, r_i)
            xof = self.init_xof(seed_e +
                int(i).to_bytes(8, byteorder='little') + b'd')
            (ct_i, k_i) = self.mmEncap_d(pp, pkl[i], r, r_i, xof)
            ctl +=  [ct_i]
            kl  +=  [k_i]
        #   endfor

        #   ct := (~ct, (~ct_i) for i in [N]
        ct  = (ctu, ctl)

        #   K := (K_i) for i in [N]
        k   = kl

        return (ct, k)


    def mmEnc_i(self, pp, r):
        """ mmEnc^i(pp; r): Shared ciphertext. """
        a   = pp
        (r_u, e_u) = r

        #   c := A * r + e_u
        r   = vec_ntt(r_u)
        c   = [ poly_zero() for _ in range(self.m)]
        for i in range(self.m):
            for j in range(self.n):
                c[i] = poly_add( c[i], mul_ntt(a[i][j], r[j]) )
        c   = vec_intt(c)
        c   = vec_add(c, e_u)

        #   u := [c mod q] (2^du)
        u   = [ poly_compress(c[i], self.du) for i in range(self.m) ]

        #   return ~ct := u
        ct  = vec_serial(u, self.du)
        return ct


    def _fast_enc(self, c):
        """ Non-randomized enc (replaces dbl, cross, compress, serial.) """
        ct  = bytearray(self.d // 8)
        k   = bytearray(self.d // 8)
        for i in range(self.d // 8):
            for j in range(8):
                x = (4 * c[8 * i + j]) // self.q
                ct[i] |= (x & 1) << j
                k[i]  |= (((x + 1) & 3) // 2) << j
        return ct, k


    def mmEncap_d(self, pp, pk_i, r, r_i, xof):
        """ mmEncap^d(pp, pk_i; r, r_i): Individual ciphertext. """

        (r_u, e_u) = r
        y   = r_i

        #   b'_i := t_i
        b,_ = vec_deserial(pk_i, self.logq, self.n)

        #   c_i := < b'_i, r > + y_i
        r   = vec_ntt(r_u)
        c   = poly_zero()
        for i in range(self.n):
            c = poly_add( c, mul_ntt(b[i], r[i]) )
        c   = intt(c)
        c   = poly_add(c, y)

        #   invoke equivalent, non-randomized encoding function
        #return self._fast_enc(c)

        #   mmEncap^d continues like this:

        #   c~_i <- dbl(c_i)
        c   = poly_dbl(c, xof)

        #   u_i <- cross( c~_i mod 2q )(2)
        u   = poly_cross(c, 2 * self.q)

        #   mu_i <- [c~_i mod 2q](2)
        mu  = poly_compress(c, 1, 2 * self.q)

        #   return (ct~_i := u_i, K_i := mu_i)
        ct  = poly_serial(u, 1)
        k   = poly_serial(mu, 1)
        return ct, k


    #   --- mmKyber:    "Short message" system

    def mmEnc(self, pp, pkl, ml, seed_e):
        """ mmEnc(pp, (pk_i), (m_i) for i in [N]): Encrypt to N recipients."""

        #   r := (r, e_u) <= D^n_sigma0 x D^M_sigma0
        xof = self.init_xof(seed_e + b'R')
        r_u = [ poly_gauss(xof, self.sigma0, self.d) for _ in range(self.n) ]
        xof = self.init_xof(seed_e + b'e')
        e_u = [ poly_gauss(xof, self.sigma0, self.d) for _ in range(self.m) ]
        r   = (r_u, e_u)

        #   ^ct <- mmEnc^i(pp; r)
        ctu = self.mmEnc_i(pp, r)

        ctl = []
        kl  = []

        #   for i in [N]
        for i in range(len(pkl)):
            #   r_i := y_i <- D_sigma1
            xof = self.init_xof(seed_e +
                int(i).to_bytes(8, byteorder='little') + b'r')
            r_i = poly_gauss(xof, self.sigma1, self.d)

            #   ~ct_i <- mmEnc^d(pp. pk_i, m_i; r, r_i)
            ct_i = self.mmEnc_d(pp, pkl[i], ml[i], r, r_i)
            ctl +=  [ct_i]
        #   endfor

        #   return ct := (~ct, (~ct_i) for i in [N]
        ct  = (ctu, ctl)
        return ct


    def mmEnc_d(self, pp, pk_i, m_i, r, r_i ):
        """ mmEnc^d(pp, pk_i, m_i; r, r_i). """
        (r_u, e_u) = r
        y   = r_i

        #   b'_i := [t_i] (q)
        b,_ = vec_deserial(pk_i, self.logq, self.n)

        #   c_i := < b'_i, r > + y_i
        r   = vec_ntt(r_u)
        c   = poly_zero()
        for i in range(self.n):
            c = poly_add( c, mul_ntt(b[i], r[i]) )
        c   = intt(c)
        c   = poly_add(c, y)

        #   .. + [q/2] * m_i
        z,_ = poly_deserial(m_i, 1)
        c   = poly_add(c, poly_decompress(z, 1))

        #   u_i := [c_i_i mod q]_2^dv
        u   = poly_compress(c, self.dv)

        #   ct_i := u_i
        ct  = poly_serial(u, self.dv)

        return ct


    def mmDec(self, pp, sk, ct):
        """ mmDec(pp, sk, ct): Decrypt a message. """

        #   deserialize
        s   = [ self.poly_nu(sk[ i*self.sz_nu : (i+1)*self.sz_nu ])
                    for i in range(self.m) ]
        (u, ui) = ct
        u,_ = vec_deserial(u, self.du, self.m)
        v,_ = poly_deserial(ui, self.dv)

        #   u' := [u mod 2^{d_v}]_2^{d_u}
        v   = poly_decompress(v, self.dv, 2**self.du)

        #   m := [ u' - <u,s> mod 2^d_u ]_2
        m   = poly_zero()
        for i in range(self.m):
            poly_mul1_add(m, u[i], s[i])
        m   =   [ (v[i] - m[i]) % 2**self.du  for i in range(256) ]
        m   = poly_compress( m, 1, 2**self.du )
        m   = poly_serial( m, 1 )
        return  m


#   paramter sets for various security levels

mmKyber_128 = mmKyber(  sec=128, N=2**10, q=MM_Q, d=256, m=4, n=4,
                        nu_bar=3, du=14, dv=2, sigma0=15.90, sigma1=368459.34)

mmKyber_192 = mmKyber(  sec=192, N=2**10, q=MM_Q, d=256, m=7, n=7,
                        nu_bar=2, du=15, dv=2, sigma0=15.90, sigma1=488797.36)

mmKyber_256 = mmKyber(  sec=256, N=2**10, q=MM_Q, d=256, m=9, n=9,
                        nu_bar=2, du=16, dv=2, sigma0=15.90, sigma1=554941.07)

