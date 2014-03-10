#! /usr/bin/env python
# Kristen K. Dang for Sage Bionetworks
# March 5, 2014

import os, subprocess, synapseclient

# Get metadata from Synapse and create a dictionary matching RNAseq data IDs to DNA IDs.
syn = synapseclient.Synapse()
syn.login()
metadataEntity = syn.get('syn2299154')
metaDict = dict()
print 'path for metafile: %s' % metadataEntity.path
with open(metadataEntity.path) as metadata:
	line_index = 1
	for line in metadata:
		if line_index == 1:
			line_index += 1
		else:
			entries = line.strip().split(',')
			print '%s %s' % (entries[47].strip('"'), entries[70])
			metaDict[entries[70].strip('"')] = entries[47].strip('"')
metadata.closed

# Iterate through VCFs and run vcf-compare on each of them vs their corresponding DNA-chip VCF data
TabixBgzips = os.listdir('/hpc/users/xdangk01/scratch/variant_concordance')
for file in TabixBgzips:
	if file.endswith('.raw.vcf.gz'):
		RNAid = file.split('.')[0]
		DNAid = metaDict[RNAid]
		outFileName = RNAid+'.compare_stats'
		
		cmd = ''.join(['vcf-compare -g -m ', DNAid, ':', RNAid, ' chip_geno_noDI.vcf.gz ', file, ' > ', outFileName])
		print '%s' % cmd
#		subprocess.call(cmd, shell = True)