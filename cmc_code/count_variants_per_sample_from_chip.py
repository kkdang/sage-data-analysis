#! /usr/bin/env python
# KKD for Sage Bionetworks
# Feb. 28, 2014

import sys, subprocess

with open(sys.argv[1], 'r') as samplesFile:
	for line in samplesFile:
		subprocess.call(' '.join(['vcf-subset -c', line.strip(), '-e /hpc/users/xdangk01/scratch/variant_concordance/chip_geno_noDI.vcf | wc -l']), shell = True)
samplesFile.closed

