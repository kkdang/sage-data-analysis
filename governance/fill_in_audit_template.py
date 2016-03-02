#! /usr/bin/env python
# Kristen K Dang for Sage Bionetworks
# August 15, 2014

import os, argparse, re, subprocess, sys
import synapseclient
import parser as pr

parser = argparse.ArgumentParser(description='Populates audit template with current audit data.')
parser.add_argument('--dateCol', required=True, help='Label of appropriate column in audit data table, e.g. AUDIT_2015_05')
parser.add_argument('--tableId', required=False, help='Synapse ID of table to process, if not default', default='syn5642510')
args = parser.parse_args()


syn = synapseclient.login()


# Read template into memory
#wiki = syn.getWiki('syn2482323', subpageId='64950') real template
#wiki = syn.getWiki('syn2482323', subpageId='67264') # test template
#template = wiki.markdown

template_entity = syn.get(entity='syn2625522', version = 9)
with open(template_entity.path, 'r') as template_file:
	template = template_file.readlines()
template_file.closed


# Get data from Synapse table
auditDataDict = dict()
results = syn.tableQuery(''.join(["select VAR_NAME,DESCRIPTION,", args.dateCol," from ", args.tableId]))
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


# When finished, print
#print '%s' % template
for line in template:
	print '%s' % line
	
	
#wiki = syn.getWiki('syn2482323', subpageId='64950') real template
# wiki = syn.getWiki('syn2482323', subpageId='67264') # test template
# wiki.markdown = 
