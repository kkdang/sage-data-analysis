#! /usr/bin/env python
# Kristen K Dang for Sage Bionetworks
# August 15, 2014

import os, argparse, re, subprocess
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
			newTemplate = populateFields(modTemplate=modTemplate,field=name,val=val)
			modTemplate = newTemplate
	users.closed
	return(modTemplate)


def populateFields(modTemplate,field,val):
	varName = ''.join(['\$',field])
	newTemplate = list()
	for line in modTemplate:
# 				if re.search(pattern=str(varName), string=str(line)):
# 					print 'match %s: %s' % (varName, line)	
		newTemplate.append(re.sub(pattern=str(varName), repl=val, string=str(line), count=0))
# 			if not re.match(pattern=str(varName), string=str(modTemplate)):
# 				print "Can't find pattern"
# 			newTemplate = re.sub(pattern=varName, repl=val, string=modTemplate, count=0)	
	return(newTemplate)


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

active_nonsage_entity = syn.get(entity=userFiles['active_non_sage_users.csv'])
ACTIVE_NONSAGE = str(int(subprocess.check_output(' '.join(['wc -l', active_nonsage_entity.path]), shell=True).split()[0])-1)
newTempl = populateFields(field='ACTIVE_NONSAGE', val=ACTIVE_NONSAGE, modTemplate=template)
template = newTempl


active_sage_entity = syn.get(entity=userFiles['active_sage_users.csv'])
ACTIVE_SAGE = str(int(subprocess.check_output(' '.join(['wc -l', active_sage_entity.path]), shell=True).split()[0])-1)
newTempl = populateFields(field='ACTIVE_SAGE', val=ACTIVE_SAGE, modTemplate=template)
template = newTempl



# Work on the PROJECTS-FILES sub-folder
projectsFiles = dict()
qr = syn.chunkedQuery(''.join(['SELECT id, name FROM entity WHERE parentId=="', existingFolders['projects_files'], '"']))
for result in qr:
	projectsFiles[result['entity.name']] = result['entity.id']
projectFilesNeeded = ['file_counts_snapshot_20141101.txt', 'file_counts_current_audit_period.txt', 'project_counts_snapshot_20141101.txt']
for dataset in projectFilesNeeded:
	newTempl = getData(filesDict=projectsFiles, fileName=dataset, modTemplate=template)
	template = newTempl



public_downloads_entity = syn.get(entity=projectsFiles['top_public_file_downloads_current_audit_period.csv'])
AUDIT_PUBLIC_DOWNLOADS = str(int(subprocess.check_output(' '.join(['wc -l', public_downloads_entity.path]), shell=True).split()[0])-1)
newTempl = populateFields(field='AUDIT_PUBLIC_DOWNLOADS', val=AUDIT_PUBLIC_DOWNLOADS, modTemplate=template)
template = newTempl

AUDIT_RESTRICTED_DOWNLOADS = subprocess.check_output(' '.join(['cut -f5 -d,', public_downloads_entity.path, '| awk \'$1 ~ /1/ {++x} END{print x}\'']), shell=True).strip()
newTempl = populateFields(field='AUDIT_RESTRICTED_DOWNLOADS', val=AUDIT_RESTRICTED_DOWNLOADS, modTemplate=template)
template = newTempl


AUDIT_CONTROLLED_DOWNLOADS = subprocess.check_output(' '.join(['cut -f4 -d,', public_downloads_entity.path, '| awk \'$1 ~ /1/ {++x} END{print x}\'']), shell=True).strip()
newTempl = populateFields(field='AUDIT_CONTROLLED_DOWNLOADS', val=AUDIT_CONTROLLED_DOWNLOADS, modTemplate=template)
template = newTempl

public_uploads_nonsage_entity = syn.get(entity=projectsFiles['public_file_uploads_non_sage_current_audit_period.csv'])
AUDIT_PUBLIC_NONSAGE_UPLOADS = str(int(subprocess.check_output(' '.join(['wc -l', public_uploads_nonsage_entity.path]), shell=True).split()[0])-1)
newTempl = populateFields(field='AUDIT_PUBLIC_NONSAGE_UPLOADS', val=AUDIT_PUBLIC_NONSAGE_UPLOADS, modTemplate=template)
template = newTempl


public_uploads_sage_entity = syn.get(entity=projectsFiles['public_file_uploads_sage_current_audit_period.csv'])
AUDIT_PUBLIC_SAGE_UPLOADS = str(int(subprocess.check_output(' '.join(['wc -l', public_uploads_sage_entity.path]), shell=True).split()[0])-1)
newTempl = populateFields(field='AUDIT_PUBLIC_SAGE_UPLOADS', val=AUDIT_PUBLIC_SAGE_UPLOADS, modTemplate=template)
template = newTempl


# When finished, print
#print '%s' % template
for line in template:
	print '%s' % line
	
	
#wiki = syn.getWiki('syn2482323', subpageId='64950') real template
# wiki = syn.getWiki('syn2482323', subpageId='67264') # test template
# wiki.markdown = 
