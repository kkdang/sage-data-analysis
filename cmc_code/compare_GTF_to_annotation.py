#! /usr/bin/env python
# KKD for Sage Bionetworks
# Feb. 11, 2014


import synapseclient, os, subprocess, pybedtools

syn = synapseclient.login()

mergedGTF_entity = syn.get(entity='syn2338565', downloadLocation = os.getcwd())
#ENSGv70_entity = syn.get(entity='syn2215531', downloadLocation = os.getcwd())

mergedGTF = pybedtools.BedTool(mergedGTF_entity.path)

#ENSGv70_exons = pybedtools.BedTool('/work/DAT_107__common_mind/Data/Homo_sapiens.GRCh37.70.processed.gtf').filter(lambda x: x[2] == 'exon').saveas('ENSGv70_exons.gff')
ENSGv70_exons = pybedtools.BedTool('ENSGv70_exons.gff')
    
def countNovelExons(fraction):
	intersectFile = ''.join(['novelExons_f', str(fraction), '.gff'])
	novelExons = mergedGTF.intersect(ENSGv70_exons, v=True, f=fraction).saveas(intersectFile)
	totalNovelCount = subprocess.call(' '.join(['awk \'{print $1"-"$4"-"$5}\'', intersectFile, '| sort | uniq | wc -l']), shell = True)
	onlyJNovelCount = subprocess.call(' '.join(['grep \'class_code "j"\'', intersectFile, '| awk \'{print $1"-"$4"-"$5}\' | sort | uniq | wc -l']), shell = True)
	print 'Required fractional match: %f; total: %d, total-J: %d' % (fraction, totalNovelCount, onlyJNovelCount)
	

countNovelExons(fraction=0.8)
