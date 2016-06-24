#! /usr/bin/env python
# KKD for Sage Bionetworks
# Jun 22, 2016

import os, argparse



parser = argparse.ArgumentParser(description='Creates input data file for metatest.')
parser.add_argument('--isoformDir', required=True, help='Directory containing FPKM files to process.')
parser.add_argument('--debug', required=False, help='Print all commands to output file even if not executed.', action='store_true')
args = parser.parse_args()


allIsoforms = dict()

allFPKMs = os.listdir(args.isoformDir)
for item in allFPKMs:
	with open(os.path.join(args.isoformDir,item), 'r') as fpkm:
		sample = item.split('.')[0]
#		print sample
		for line in fpkm:
			if line.startswith('ENST'):
				vals = line.strip().split()
				if vals[0] in allIsoforms:
					allIsoforms[vals[0]][sample] = (vals[9],vals[10],vals[11],vals[12])
				else:
					allIsoforms[vals[0]] = {sample: (vals[9],vals[10],vals[11],vals[12])}
	fpkm.close()


for transcript in allIsoforms:
	for item in allFPKMs:
		sample = item.split('.')[0]		
		if sample in allIsoforms[transcript]:
#			print sample, allIsoforms[transcript][sample]
			lineout = '%s\t%s\t%s\t%s\t%s\t%s' % (transcript, allIsoforms[transcript][sample][0], allIsoforms[transcript][sample][1], allIsoforms[transcript][sample][2], allIsoforms[transcript][sample][3], sample)
		else:
			lineout = '%s\tNA\tNA\tNA\tNA\t%s' % (transcript, sample)
		print lineout