#	make_tables.py
#	===	parse data and create contents of tables

import sys

args = ' '.join(sys.argv[1:])

alg_l = set()
var_l = set()
data = {}

#	parse benchmark data from stdin

for line in sys.stdin:

	if line == '' or line[0] == '#':
		continue

	v = line.split()
	if len(v) < 6:
		continue

	alg = v[0]
	alg_l |= { alg }
	var	= v[1]
	var_l |= { var }

	if v[2] == 'N=':
		n = int(v[3])
	else:
		continue

	if v[4] == 'len=':
		val = int(v[5])				# bytes
	elif v[4] == 'cyc=':
		val = int(v[5])				#	cycles

		if var != 'mmSetup()':
			val /= n
		
		if len(v) >= 8 and v[6] == 'sec=':
			sec = float(v[7])
		else:
			sec = None
	else:
		continue

	if alg not in data:
		data[alg] = {}
	if var not in data[alg]:
		data[alg][var] = {}
	if n in data[alg][var]:
		print(f'double item {alg} {var} {n}')
	data[alg][var][n] = val


#	columns

kyber_l	=	[ f'Kyber{4*sec}' for sec in [ 128, 192, 256 ] ]
mmpke_l	=	[ f'mmKyber-PKE-{sec}' for sec in [ 128, 192, 256 ] ]
mmkem_l	=	[ f'mmKyber-KEM-{sec}' for sec in [ 128, 192, 256 ] ]

#	for non-batch list

non_batch 	= [	( 'mmSetup()',	mmkem_l ),
				( 'mmKGen()',	mmkem_l ),
				( 'mmDec()',	mmpke_l ),
				( 'mmDecap()',	mmkem_l ),
				( 'ML-KEM.KeyGen()', kyber_l ),
				( 'ML-KEM.Decaps()', kyber_l ),
				( 'K-PKE.Decrypt()', kyber_l )	]

for	var in non_batch:
	s	= f'\t{{\\sf {var[0]:16}}} '
	for alg in var[1]:
		s += f' & {data[alg][var[0]][1024]:10,.0f}'
	s	+= ' \\\\'
	print(s)

s = '{\\bf Algorithm}'
for nn in range(0, 11):
	s += f' & $N=2^{{{nn}}}$'
	if nn == 4:
		s += '\n'
s += ' \\\\'
print(s)

for sec in [128, 192, 256]:
	alg	= [ f'mmKyber-PKE-{sec}', f'mmKyber-KEM-{sec}', f'Kyber{4*sec}' ]
	lab	= [ f'mmKyber-PKE-{sec}', f'mmKyber-KEM-{sec}', f'ML-KEM-{4*sec}' ]
	var = [ 'mmEnc()', 'mmEncap()', 'K-PKE.Encrypt()' ]
	print('\\midrule')
	for i in range(3):
		s	= f'{{\\sf {lab[i]:15} }}'
		for nn in range(0, 11):
			s += f' & {data[alg[i]][var[i]][2**nn]:7,.0f}'
			if nn == 4:
				s += f'\n{'':12}'
		s += ' \\\\'
		print(s)
	
