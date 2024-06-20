#!/usr/bin/python3.6

import sys
from optparse import OptionParser

parser = OptionParser(usage='Usage: '+sys.argv[0]+'/path/to/id.org.tab, the file created by parseGBFF.py', description='Parses output of parseGBFF.py')
parser.add_option("-d", "--debug", dest="doDebug", action='store_true', default=False, help='debug')
parser.add_option("-1", "--header", dest="outHeader", action='store_true', default=False, help='Add header 1st row')

(options, args) = parser.parse_args(args=None, values=None)
if options.doDebug: print ('-----Options------'); print(args); print(options)

#print(options)
#print(args)

if len(args) != 1:
	parser.print_help()
	sys.exit()

import pandas
import datetime
import json
from dateutil.parser import parse
#from zipfile import ZipFile


def CustomParser(data):
	j1 = json.loads(data)
	return j1

df = pandas.read_csv(args[0],sep='\t',dtype={1:str},converters={4:CustomParser},header=None)
epoch = parse('1000-01-01')
for i,row in df.iterrows():
	if options.doDebug: print (i,row)
	df.loc[i,'pathogen_tax_id']=row[1]
	df.loc[i,'pathogen']=row[4]['organism'] if 'organism' in row[4] else None
	df.loc[i,'collection_date']=parse(row[4]['collection_date'],fuzzy=True, default=epoch) if 'collection_date' in row[4] else None
	df.loc[i,'country']=row[4]['country'] if 'country' in row[4] else None
	df.loc[i,'host']=row[4]['host'] if 'host' in row[4] else None
	df.loc[i,'serotype']=row[4]['serotype'] if 'serotype' in row[4] else None
	if 'strain' in row[4]:
		df.loc[i,'isolate']=row[4]['strain']
	elif 'isolate' in row[4]:
		df.loc[i,'isolate']=row[4]['isolate']

#	df.loc[i,'strain']=row[4]['strain'] if 'strain' in row[4] else None
df_ = df[~df['collection_date'].isnull() & ~df['country'].isnull() & ~df['host'].isnull()]
#print(df)
#print(df_)
df_.to_csv(sys.stdout, columns=['pathogen_tax_id','pathogen','collection_date','country','host','serotype','isolate'], index=False, header=options.outHeader)

# FF=open('filtered.fasta','w')
# for filename in ('fasta/DNA/'+df_[0]).tolist():
# 	f=open(filename,'r')
# 	FF.write(f.read())
# 	f.close()
# FF.close()
#	f=open('fasta/DNA/'+row[0])
