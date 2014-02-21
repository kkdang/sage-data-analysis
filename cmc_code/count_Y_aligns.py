#! /usr/bin/env python
# KKD for Sage Bionetworks
# Feb. 5, 2014

import os, sys, subprocess

tophatDir = '/projects/CommonMind/data/FROM_CORE/Production/fastQ/BySample/tophat209.2.2'

bamDirList = os.listdir(tophatDir)
for bamDir in bamDirList:
	filesList = os.listdir(os.path.join(tophatDir, bamDir))
	for file in filesList:
		if file.endswith('.sort.bam'):
			cmd = ' '.join(['samtools view -c', os.path.join(tophatDir, bamDir, file), '\'chrY\''])
			print '%s' % cmd
#			subprocess.call(cmd, shell = True)