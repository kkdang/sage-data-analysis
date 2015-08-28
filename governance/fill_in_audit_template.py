#! /usr/bin/env python
# Kristen K Dang for Sage Bionetworks
# August 15, 2014

import os, argparse, re, subprocess, sys
import synapseclient
import parser as pr

parser = argparse.ArgumentParser(description='Populates audit template with current audit data.')
parser.add_argument('--currentFolderId', dest='fid', required=True, help='Synapse ID of folder containing audit data files to process.')
parser.add_argument('--dateCol', required=True, help='Label of appropriate column in audit data table, e.g. AUDIT_2015_05')
args = parser.parse_args()


syn = synapseclient.login()


# Read template into memory
#wiki = syn.getWiki('syn2482323', subpageId='64950') real template
#wiki = syn.getWiki('syn2482323', subpageId='67264') # test template
#template = wiki.markdown

template_entity = syn.get(entity='syn2625522', version = 8)
with open(template_entity.path, 'r') as template_file:
	template = template_file.readlines()
template_file.closed


# Get data from Synapse table
auditDataDict = dict()
results = syn.tableQuery(''.join(["select VAR_NAME,DESCRIPTION,", args.dateCol," from syn4892864"]))
for result in results:
	auditDataDict[result[2]] = float(result[4])
   	print '%s = %s' % (result[2], result[4])



# Read math in, calculate and put in dictionary
with open('/Users/kristen/Computing/Audit/audit_formulas.txt', 'r') as formulas_file:
	for line in formulas_file:
		(valName,formula) = line.split('\t')
		compFormula = pr.expr(formula).compile()
		auditDataDict[valName] = eval(compFormula)
		print '%s = %s' % (valName, eval(compFormula))
formulas_file.closed


def populateFields(modTemplate,field,val):
	varName = ''.join(['\$',field])
	newTemplate = list()
	for line in modTemplate:
		newTemplate.append(re.sub(pattern=str(varName), repl=val, string=str(line), count=0))
	return(newTemplate)


# Replace variables in text with values in dictionary (first cast values to string)
modTemplate = template
for field in auditDataDict:
	newTemplate = populateFields(modTemplate=modTemplate,field=field,val=str(auditDataDict[field]))
	modTemplate = newTemplate
template = modTemplate

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
