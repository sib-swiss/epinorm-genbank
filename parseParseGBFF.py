#!/usr/bin/python3.6

import sys
import os
from optparse import OptionParser

import pandas
import datetime
import json
from dateutil.parser import parse
#from zipfile import ZipFile

parser = OptionParser(usage='Usage: '+sys.argv[0]+'/path/to/id.org.tab [/path/to/NCBItaxNames.csv]', description='Parses id.org.tab, the file created by parseGBFF.py; Optional file NCBItaxNames.csv contains translation of host names to NCBI taxids')
parser.add_option("-1", "--header", dest="outHeader", action='store_true', default=False, help='add 1st row with column names')
parser.add_option("-i", "--ignore", dest="imed", action='store_true', default=False, help='ignore absence of any of 3 outbreak variables: geo, date, host species')
parser.add_option("-d", "--debug", dest="doDebug", action='store_true', default=False, help='debug')

(options, args) = parser.parse_args(args=None, values=None)
if options.doDebug: print ('-----Options------'); print(args); print(options)

if len(args) == 0:
	parser.print_help()
	sys.exit()

# to parse JSON data
def CustomParser(data):
	j1 = json.loads(data)
	return j1
#
# read main data
df = pandas.read_csv(args[0],sep='\t',dtype={1:str},converters={4:CustomParser},header=None)

epoch = parse('1000-01-01')
for i,row in df.iterrows():
	if options.doDebug: print (i,row)
	df.loc[i,'Pathogen NCBI taxonomy ID']=row[1]
	df.loc[i,'Pathogen species']=row[4]['organism'] if 'organism' in row[4] else None
	try:
		df.loc[i,'Date observed']=parse(row[4]['collection_date'], default=epoch, fuzzy=True) if 'collection_date' in row[4] else None
	except ValueError:
		sys.stderr.write('Cannot parse date: '+row[4]['collection_date']+' in record: '+ str(i)+", the record skipped\n" )
		continue
	df.loc[i,'Geo text original']=row[4]['country'] if 'country' in row[4] else None
	df.loc[i,'host']=row[4]['host'].lower() if 'host' in row[4] else None
	df.loc[i,'Host species Latin name']=row[4]['host'] if 'host' in row[4] else None
	df.loc[i,'Pathogen serotype']=row[4]['serotype'] if 'serotype' in row[4] else None
	if 'strain' in row[4]:
		df.loc[i,'Pathogen isolate or strain']=row[4]['strain']
	elif 'isolate' in row[4]:
		df.loc[i,'Pathogen isolate or strain']=row[4]['isolate']

# if present, read taxonomy data and join it to main data using lowercased sci name
if len(args)==2 and args[1] is not None and os.path.isfile(args[1]):
	tx = pandas.read_csv(args[1],dtype={'tax_id':str})
	df = pandas.merge(df,tx,how='left',left_on='host',right_on='synonym_name')
	df['Host species NCBI taxonomy ID'] = df['tax_id']
else:
	sys.stderr.write('No data to evaluate host species NCBI tax ids\n')
	df['Host species NCBI taxonomy ID'] = None

# filter out records having no either of geo, date or host species data
if not options.imed:
	df = df[~df['Date observed'].isnull() & ~df['Geo text original'].isnull() & ~df['Host species Latin name'].isnull()]

# print out selected columns
df.to_csv(sys.stdout, columns=['Pathogen NCBI taxonomy ID','Pathogen species','Pathogen serotype','Pathogen isolate or strain','Host species Latin name','Host species NCBI taxonomy ID','Date observed','Geo text original'], index=False, header=options.outHeader)

### optionally dump corresponding fasta files.
# FF=open('filtered.fasta','w')
# for filename in ('fasta/DNA/'+df_[0]).tolist():
# 	f=open(filename,'r')
# 	FF.write(f.read())
# 	f.close()
# FF.close()
#	f=open('fasta/DNA/'+row[0])



