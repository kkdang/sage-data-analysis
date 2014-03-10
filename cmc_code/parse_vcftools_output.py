#! /usr/bin/env python
# Kristen K. Dang for Sage Bionetworks
# March 6, 2014

import os

dataPath = '/hpc/users/xdangk01/scratch/variant_concordance'
statsFiles = subprocess.check_output(' '.join(['find', dataPath, '-name \'*.compare_stats\'']))
for file in statsFiles:
	vals = list()
	printLine = os.path.basename(file.strip())
	for line in subprocess.check_output(' '.join(['grep ^[SV]N', file, '| cut -f 2-'])):
		items = line.strip().split()
		if line.startswith('[0-9]'):
			printLine+'\t'+line[0]
#			vals.append(line[0])
		elif line.startswith('')
			printLine+'\t'+line[1]
#			vals.append(line[1])
#	for val in vals:
#		printLine+'\t'+val
	print '%s' % printLine