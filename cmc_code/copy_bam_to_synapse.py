import synapseclient, os
from synapseclient import File

syn = synapseclient.login()

# Make dictionary of filenames to load
releasedFiles = dict()
with open('/hpc/users/xdangk01/CMC_released_RNAseq') as released:
	for line in released:
		releasedFiles[line.strip()] = 1

# Make dictionary of files already loaded		
loadedEntities = dict()
results = syn.chunkedQuery(''.join(['SELECT id, name FROM entity WHERE parentId=="syn4645130"']))
for result in results:
	loadedEntities[result['entity.name']] = result['entity.id']

# Iterate through all files, loading those that are intended for release but are not yet loaded.
with open('/sc/orga/projects/CommonMind/data/FROM_CORE/Phase1/BAMS.fofn') as bamFiles:
	for line in bamFiles:
		if os.path.basename(line.strip()).split('.')[0] in releasedFiles and not os.path.basename(line.strip()) in loadedEntities:
			syn.store(File(path=line.strip(),parent='syn4645130',synapseStore=True))
			print 'copied %s' % line.strip()
		else:
			print '%s' % line.strip()
bamFiles.close()
