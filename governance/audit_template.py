#! /usr/bin/env python
# Kristen K Dang for Sage Bionetworks
# August 15, 2014

import os, argparse, re
import synapseclient

parser = argparse.ArgumentParser(description='Populates audit template with current audit data.')
parser.add_argument('--currentFolderId', dest='fid', required=True, help='Synapse ID of folder containing audit data files to process.')
args = parser.parse_args()


syn = synapseclient.login()


# Read template into memory
#wiki = syn.getWiki('syn2482323', subpageId='64950') real template
#wiki = syn.getWiki('syn2482323', subpageId='67264') # test template
#template = wiki.markdown

template_entity = syn.get(entity='syn2625522')
with open(template_entity.path, 'r') as template_file:
	template = template_file.readlines()
template_file.closed


def getData(filesDict,fileName,modTemplate=template):
	'''get data in file and populate the template with that data'''

	entity = syn.get(entity=filesDict[fileName])
	data = dict()
	with open(entity.path, 'r') as users:
		for line in users:
			(name, val, description) = line.split('\t')
			varName = ''.join(['\$',name])
			newTemplate = list()
			for line in modTemplate:
# 				if re.search(pattern=str(varName), string=str(line)):
# 					print 'match %s: %s' % (varName, line)	
				newTemplate.append(re.sub(pattern=str(varName), repl=val, string=str(line), count=0))
# 			if not re.match(pattern=str(varName), string=str(modTemplate)):
# 				print "Can't find pattern"
# 			newTemplate = re.sub(pattern=varName, repl=val, string=modTemplate, count=0)	
			modTemplate = newTemplate
	users.closed
	return(modTemplate)



# Get Synapse IDs of sub-folders
existingFolders = dict()
results = syn.chunkedQuery(''.join(['SELECT id, name, concreteType FROM entity WHERE parentId=="', args.fid, '"']))
for result in results:
	if result['entity.concreteType'][0] == 'org.sagebionetworks.repo.model.Folder':
		existingFolders[result['entity.name']] = result['entity.id']


# Work on the users sub-folder
userFiles = dict()
qr = syn.chunkedQuery(''.join(['SELECT id, name FROM entity WHERE parentId=="', existingFolders['users_activities'], '"']))
for result in qr:
	userFiles[result['entity.name']] = result['entity.id']

userFilesNeeded = ['user_counts.txt']

for dataset in userFilesNeeded:
	newTempl = getData(filesDict=userFiles, fileName=dataset)
	template = newTempl


# When finished, print
#print '%s' % template
for line in template:
	print '%s' % line