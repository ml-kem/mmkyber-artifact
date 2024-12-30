#   mm_256_param.py
#   binary secret

vname   = "mm_256_param"    #   variable name

deg     = 256               #   ring Rp degree d
mod     = 33550337          #   ring Rp modulus p
m,n     = 9,18              #   matrix [ A | I ]
dim     = (m,n)             #   dimensions of A in Rp^(m,n)

# partition of w : [w]
wpart   = [ list(range(0,9)), list(range(9,18)) ]

#wl2    = [ 1.2 * ((1/2)*n*deg)**0.5 ]
wl2     = [ 0, 0 ]          #   l2 norm (no)
wbin    = [ 1, 1 ]          #   binary (yes)
#rej    = [ 1 ]             #   rej. sampling

wlinf   = 1                 #   optional "helper" linf

