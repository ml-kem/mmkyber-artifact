#	table_rel.py
#	===	parse data and create contents of tables

#	some relative data..

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

nn	=	1024

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
		s += f' & {data[alg][var[0]][nn]:10,.0f}'
	s	+= ' \\\\'
	print(s)

for sec in [128, 192, 256]:
	alg	= [ f'mmKyber-KEM-{sec}', f'mmKyber-PKE-{sec}', f'Kyber{4*sec}' ]
	lab	= [ f'mmKyber-KEM-{sec}', f'mmKyber-PKE-{sec}', f'ML-KEM-{4*sec}' ]
	var = [ 'mmEncap()', 'mmEnc()', 'K-PKE.Encrypt()' ]
	ct	= 'ciphertext'
	print('\\hline')
	for i in range(3):
		s = f'{{\\sf {lab[i]:15} }}'
		s += f' & {data[alg[i]][var[i]][nn]:7,.0f}' 
		s += f' & {data[alg[i]][ct][nn]:9,.0f}'
		if i < 2:
			x = data[alg[2]][var[2]][nn] / data[alg[i]][var[i]][nn]
			s += f' & {x:5.2f}'
			x = data[alg[2]][ct][nn] / data[alg[i]][ct][nn]
			s += f' & {x:5.1f}'
		print(s)
	