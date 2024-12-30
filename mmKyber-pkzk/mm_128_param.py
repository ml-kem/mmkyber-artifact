#   mm_128_param.py
#   ternary secret

vname   = "mm_128_param"    #   variable name

deg     = 256               #   ring Rp degree d
mod     = 33550337          #   ring Rp modulus p
m,n     = 4,16              #   matrix [ A | I | -A | -I ]
dim     = (m,n)             #   dimensions of A in Rp^(m,n)

# partition of w : [w]
wpart   = [ list(range( 0, 4)), list(range( 4, 8)),
            list(range( 8,12)), list(range(12,16)) ]

#wl2    = [ 1.2 * ((2/3)*n*deg)**0.5 ]
wl2     = [ 0, 0, 0, 0 ]    #   l2 norm (no)
wbin    = [ 1, 1, 1, 1 ]    #   binary (two times)
#rej    = [ 1 ]             #   rej. sampling
wlinf   = 1                 #   optional "helper" linf

