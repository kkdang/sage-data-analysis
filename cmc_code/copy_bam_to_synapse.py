import synapseclient, os
from synapseclient import File

syn = synapseclient.login()

releasedFiles = dict()
with open('/hpc/users/xdangk01/CMC_released_RNAseq') as released:
	for line in released:
		releasedFiles[line.strip()] = 1
		
loadedEntities = dict()
results = syn.chunkedQuery(''.join(['SELECT id, name FROM entity WHERE parentId=="syn4922873"']))
for result in results:
	loadedEntities[result['entity.name']] = result['entity.id']


with open('/sc/orga/projects/CommonMind/data/FROM_CORE/Phase1/UNMAPPED_BAMS.fofn') as bamFiles:
	for line in bamFiles:
		sampleName = line.strip().split('/')[11].lstrip('tophat_')
		sampleFile = '.'.join([sampleName, 'unmapped.bam'])
			
		if sampleName in releasedFiles and not sampleFile in loadedEntities:
			syn.store(File(path=line.strip(),parent='syn4922873',synapseStore=True,name=sampleFile))
			print 'copied %s' % sampleFile
		else:
			print '%s' % sampleFile
bamFiles.close()
