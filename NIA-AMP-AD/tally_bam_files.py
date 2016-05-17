#! /usr/bin/env python

import synapseclient

syn = synapseclient.login()


seqFolders = {'taud35':'syn5759582', 'tau-mouse-app':['syn4486837', 'syn4486995'], 'mssm-mouse':'syn5053456', 'mssm-drosophila':['syn5012222','syn5012145'], 'microglia':'syn5053443', 'broad-ipsc':'syn4228582', 'broad-mdmi':'syn4228560', 'rosmap':'syn4164376', 'mayo-wes':['syn5012301','syn5519740'], 'mssm':['syn4883087', 'syn5480006', 'syn4055270'], 'mayo-pilot-psp':['syn5584594', 'syn4518661'], 'mayo-pilot-ad':['syn5580964','syn4518517'], 'mayo-tcx':'syn4894912', 'mayo-cbe':'syn5049322'}


def getFileSizes(folderId):
	results = syn.chunkedQuery(''.join(['select id from file where parentId=="', folderId, '"']))
	sumSize = 0
	for result in results:
#		print result['file.id']
		ent = syn.get(result['file.id'], downloadFile = False)
		sumSize += int(ent.fileSize)
	return(sumSize)


total = 0
for key,val in seqFolders.iteritems():
	print key
	if isinstance(val, list):
		for element in val:
			total += getFileSizes(element)			
	else:
		total += getFileSizes(val)
		
print total/1e9