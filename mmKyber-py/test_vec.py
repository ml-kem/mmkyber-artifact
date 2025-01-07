#   test_vec.py
#   === Test Vector Generation

from mmKyber import mmKyber_128, mmKyber_192, mmKyber_256

#   for creating "header files" for LaZer experiments
#from test_lazer import lazer_h

if (__name__ == "__main__"):

    def dbg_sum(s, lab=''):
        """ Simple polynomial checksum. """
        x = 1
        for b in s:
            x = (x * 0x103 + b) % (2**32)
        print(f'{lab}[{len(s)}] chk {x:08x}')
        return x

    def testvec(iut, is_kem, N):

        seed_a = b'0123456789abcdef'
        seed_k = b'000102030405060708090a0b0c0d0e0f'
        seed_e = b'00112233445566778899aabbccddeeff'

        if is_kem:
            alg = f'mmKyber-KEM-{iut.sec}'
        else:
            alg = f'mmKyber-PKE-{iut.sec}'

        print(f'=== chk {alg} N={N}')

        pp  = iut.mmSetup(seed_a)

        pkl = []
        skl = []
        for i in range(N):
            seed_k      = i.to_bytes(8, byteorder='little') + seed_k[8:32]
            pk_i, sk_i  = iut.mmKGen(pp, seed_k)
            pkl += [pk_i]
            skl += [sk_i]
            dbg_sum(sk_i, 'sk')
            dbg_sum(pk_i, 'pk')

            #   write the lazer header file
            #lazer_h(iut, f'{alg}-data-{i}.h', pp, pk_i, sk_i)

        if is_kem:
            #   encapsulation
            ct, kl  = iut.mmEncap(pp, pkl, seed_e)
        else:
            #   encryption
            ml  =   []
            for i in range(N):
                buf = seed_e + i.to_bytes(8, byteorder='little') + b'm'
                ml += [ iut.xof_a.new(buf).read(32) ]
            ct  =   iut.mmEnc(pp, pkl, ml, seed_e)
        (ctu, ctl) = ct
        dbg_sum(ctu, 'ct_u')

        ct = ctu + b''.join(ctl)
        dbg_sum(ct, 'ct')

        for i in range(N):
            if is_kem:
                k_i = iut.mmDecap(pp, skl[i], (ctu, ctl[i]))
                print(f'decaps #{i}: {k_i == kl[i]}')
                dbg_sum(kl[i], 'k_i')
            else:
                m_i = iut.mmDec(pp, skl[i], (ctu, ctl[i]))
                print(f'dec #{i}: {m_i == ml[i]}')
                dbg_sum(m_i, 'm_i')
            dbg_sum(ctl[i], 'ct_i')


    #   test targets
    testvec( iut=mmKyber_128, is_kem=True, N=5 )
    testvec( iut=mmKyber_128, is_kem=False, N=5 )
    testvec( iut=mmKyber_192, is_kem=True, N=5 )
    testvec( iut=mmKyber_192, is_kem=False, N=5 )
    testvec( iut=mmKyber_256, is_kem=True, N=5 )
    testvec( iut=mmKyber_256, is_kem=False, N=5 )

