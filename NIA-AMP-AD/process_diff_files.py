#! /usr/bin/env python
# KKD for Sage Bionetworks
# Apr 18, 2016


import os, sys


firstFileVals = dict()
allDiffFiles = os.listdir(sys.argv[1])
for dataFile in allDiffFiles:
	if dataFile.endswith('.diffs'):
		with open(dataFile, 'r') as diffResults:
			for line in diffResults:
				if line.startswith('< ENSG'):
					vals1 = line.lstrip('< ').split()
					firstFileVals[vals1[0]] = vals1[6]
				elif line.startswith('> ENSG'):
					vals2 = line.lstrip('> ').split()
					if vals2[0] in firstFileVals:
						diffVal = int(firstFileVals[vals2[0]]) - int(vals2[6])
						meanVals = (int(firstFileVals[vals2[0]]) + int(vals2[6])) / 2
						print '%s\t%s\t%d\t%f' % (dataFile, vals2[0], diffVal, meanVals)
						del firstFileVals[vals2[0]]
					else:
						print 'Could not find a match for '+vals2[0]+'\n'
