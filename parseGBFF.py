#!/usr/bin/python3 -u
import sys
#import time
import glob
import os
import re
import base64
import hashlib
import json
import gzip
from Bio import SeqIO
from Bio import Entrez
from collections import OrderedDict
from collections import defaultdict
from optparse import OptionParser

parser = OptionParser(usage='Usage: parseGBFF.py [options] /path/to/genebank/gzippedfiles', description='Parses gzipped GenBank file format (GBFF) files into a set of tab-delimited files, as well as fasta files into auto-created subdir fasta/. Requires a path to GB files as argument.')
parser.add_option("--outputdir", dest="outDir", type='string', default='.', help='path where output files are created, the current dir by default.')
parser.add_option("--fileformat", dest="srcFileFormat", type='string', default='genbank', help='source data file format, supported values "genbank" [default] or "embl"')
parser.add_option("--sourcefeatures", dest="mySourceFeatures", type='string', default='OS', help='source features to extract, in form of a string of uppercased characters (default OS), according to O=>organism, S=>strain, I=>isolate, T=>serotype, H=>host, G=>country, D=>collection_date, C=>clone')
parser.add_option("--seqids", dest="filterSeqIdFile", type='string', help='optional file with GenBank IDs to filtered against; one ID per line, with or without version.')
parser.add_option("--taxids", dest="filterTaxIdFile", type='string', help='optional file with NCBI tax IDs to filtered against; one ID per line.')
parser.add_option("--complete", dest="isComplete", action='store_true', default=False, help='flag to require presence of "complete xxx" (xxx=CDS, genome, sequence) in GB entry description')
parser.add_option("--cds", dest="mustHaveCDS", action='store_true', default=False, help='flag to require presence of CDS section(s) in GB entry table of features')
parser.add_option("--geneforcds", dest="mustHaveGeneForCDS", action='store_true', default=False, help='flag to require usable gene section to process CDS section(s) in GB entry table of features')
parser.add_option("--musttaxon", dest="mustHaveTaxon", action='store_true', default=False, help='flag to require presence of NCBI taxon and organism name in GB entry source section')
parser.add_option("--levelncbitaxid", dest="topNCBItaxon", type='string', default=None, help='check each parsed NCBI taxons is under this clade taxon, e.g. 2 (Bacteria); leave blank otherwise')
parser.add_option("--excludedups", dest="excludeAAIDdups", action='store_true', default=False, help='exclude [take the first and ignore remaining ones] duplicating genomic IDs within same organism ID')
parser.add_option("--zapclade", dest="excludeTopNCBItaxon", type='string', default=None, help='check each parsed NCBI taxons is NOT under this clade taxon, e.g. 2 (Bacteria); leave blank otherwise')
parser.add_option("--redoTaxon", dest="redoTaxon", type='string', default='', help='no value - do nothing; Y - create pseudo taxon aka organism id as <taxon>_0; DB - ask for its value from DB')
parser.add_option("--dnafasta", dest="dnaFasta", action='store_true', default=False, help='make FASTA files with DNA sequences read from GB entries')
parser.add_option("--aafasta", dest="aaFasta", action='store_true', default=False, help='make FASTA files with AA sequences read from GB entries')
parser.add_option("--fastafilename", dest="fastaFilename", type='string', default='D', help='fastafile file name if interpretation mode is G: D - substr 1:10 (md5(JSON(organism description))) [default]; A - Assembly ID; a - Assembly ID without version suffix; T - biggest NCBI tax id from record [optionally verified, see -l]. If interpretation mode is A, fasta file name is md5 of the GBFF source file name')
parser.add_option("--allseqfiles", dest="allSeqFiles", action='store_true', default=False, help='dump DNA and AA sequences in both non-redundant and redundant files, hence total 4 files')
parser.add_option("--fetchfromentrez", dest="fetchFromENTREZ", action='store_true', default=False, help='Online-fetch absent protein data [of polyprotein origin] using Entrez e-utils')
parser.add_option("--interpretemode", dest="interpreteMode", type='string', default='A', help='records interpretation mode: A - assembly file, i.e. we keep the tax id from the first record for the entire assembly file [default]; G - GeneBank file, use tax ids from each record')
parser.add_option("--gidmode", dest="gidMode", type='string', default='G', help='G (default) - include gid at first priority as local gene id, N - include gene name at first priority as local gene id')
parser.add_option("--debug", dest="doDebug", action='store_true', default=False, help='print BioPython representation of GB record into stdout')

(options, args) = parser.parse_args(args=None, values=None)
if options.doDebug: print ('-----Options------'); print(args); print(options)

#print(options)
#print(args)

if len(args) != 1:
	parser.print_help()
	sys.exit()
	#sys.exit("Usage: " +sys.argv[0]+ " /path/to/genebank/gzippedfiles /path/to/output/files [filter.file.with.ids]")
else:
	gbFilesDir = args[0]  
	if not os.path.exists(gbFilesDir):
		sys.exit(gbFilesDir + " does not exist")
	outDir = options.outDir #sys.argv[2] if len(sys.argv) == 3 else '.'
	if not os.path.exists(outDir):
		sys.exit(outDir + " does not exist")
	outFastaDNAdir = outDir+'/fasta/DNA'
	outFastaAAdir  = outDir+'/fasta/AA'
	if options.dnaFasta:
		if not os.path.exists(outFastaDNAdir):
		    os.makedirs(outFastaDNAdir)
	if options.aaFasta:
		if not os.path.exists(outFastaAAdir):
		    os.makedirs(outFastaAAdir)


fpIdsDNA=open(outDir+'/id.DNA.tab','w')
fpIdsAA =open(outDir+'/id.AA.tab','w')
fpIdsAAadd =open(outDir+'/id.AA.add.tab','w')
fpOrgIDs=open(outDir+'/id.org.tab','w')

if options.allSeqFiles: fpSeqNrDNA=open(outDir+'/seq.nr.DNA.tab','w')
else: fpSeqNrDNA = None
if options.allSeqFiles: fpSeqAllDNA=open(outDir+'/seq.all.DNA.tab','w')
else: fpSeqAllDNA = None
if options.allSeqFiles: fpSeqAA =open(outDir+'/seq.nr.AA.tab','w')
else: fpSeqAA = None
#fpFstDNA=open(outDir+'/seq.DNA.fasta','w')                             ;else: fpFstDNA = None;
#fpKW =open(outDir+'/kw.tab'  ,'w')

DNAseqHash=set(); AAseqHash=set(); SeqIdsToFilter=set(); orgIDhash = set();

noFilter=True
if options.filterSeqIdFile != None:
	noFilter = False
	SeqIdsToFilter = set(open(options.filterSeqIdFile, 'rt').read().split()) # set is MUCH faster than list
	fpAbsentIDs=open(outDir+'/id.DNA.absent.tab','w')

if options.filterTaxIdFile != None:
	TaxIdsToFilter = set(map(int,open(options.filterTaxIdFile, 'rt').read().split())) # set is MUCH faster than list

def write_into_taxon_file(a_dir,fastafilename,sequence,idstr):
	if options.dnaFasta or options.aaFasta:
		fp = open(a_dir+'/'+fastafilename, 'at')
		fp.write('>'+idstr+"\n")
		fp.write(sequence+"\n")
		fp.flush()
		fp.close()
	return True
def get_gb_gualifiers(a_list,cv):
	d = {}
	for r in a_list:
		if r.get('GBQualifier_value') is not None and r.get('GBQualifier_name') is not None:
			if str(r['GBQualifier_name']) in cv:
				d[r['GBQualifier_name']] = str(r['GBQualifier_value'])
	return (d)
polyprotein_cv = {'Region':{'region_name','note'},'mat_peptide':{'product','peptide'}}

#if options.redoTaxon is not None: allTaxons = defaultdict(set)
#if options.redoTaxon is not None: allTaxonsPerFile = defaultdict(set)

def getAllTaxonsInClade(tx):
	cursor.execute(
		"with _1 as (select leaftaxon FROM dwnld.ncbitaxtreeinverse "
		"where array_position(txpath," + tx + ") is not null) "
		"select leaftaxon as tax_id from _1 "
		"union "
		"SELECT ncbitaxnodesmerged.tax_id FROM _1 "
		"join dwnld.ncbitaxnodesmerged on "
		"_1.leaftaxon=ncbitaxnodesmerged.tax_id_new "
	)
	return(set(i[0] for i in cursor.fetchall()))


# if topNCBItaxon given, prepare beforehand a hash table of all matching taxon
# much faster vs. consulting DB for each one
if options.topNCBItaxon is not None and options.topNCBItaxon > '':
	#cursor.execute(
	#	"with _1 as (select leaftaxon FROM dwnld.ncbitaxtreeinverse "
	#	"where array_position(txpath," + options.topNCBItaxon + ") is not null) "
	#	"select leaftaxon as tax_id from _1 "
	#	"union "
	#	"SELECT ncbitaxnodesmerged.tax_id FROM _1 "
	#	"join dwnld.ncbitaxnodesmerged on "
	#	"_1.leaftaxon=ncbitaxnodesmerged.tax_id_new "
	#)
	#allTaxonsUnderTop = set(i[0] for i in cursor.fetchall())
	allTaxonsUnderTop = getAllTaxonsInClade(options.topNCBItaxon)
	options.mustHaveTaxon = True
	#if options.doDebug: print ('-----must NCBI taxons------\n'); print(allTaxonsUnderTop)

if options.excludeTopNCBItaxon is not None and options.excludeTopNCBItaxon > '':
	allExcludedTaxonsUnderTop = getAllTaxonsInClade(options.excludeTopNCBItaxon)


#my_source_features = ['organism','strain','isolate','serotype','host','segment']
#my_source_features = ['organism','strain','isolate','serotype','host','country','collection_date']#,'clone']
all_source_features = OrderedDict({'O':'organism','S':'strain','I':'isolate','T':'serotype','H':'host','G':'country','D':'collection_date','C':'clone'}) # OSITHGDC
C = set()
for c in options.mySourceFeatures : C.add(c)
a = dict()
a = OrderedDict({key:value for key,value in all_source_features.items() if key in C})
my_source_features = list(a.values())
if options.doDebug: print ('-----Record------'); print(my_source_features)
flagComplete = re.compile('complete\s+(\w+\s+){0,2}(genome|cds|sequence|segment|rna|dna)', re.IGNORECASE)
flagWGS = re.compile('wgs|(whole genome shotgun sequence)', re.IGNORECASE)
flagNPatch = re.compile('NOVEL_PATCH')
flagAltLoc = re.compile('ALTERNATE_LOCUS')

seqfileOrgIDs = dict()

if options.srcFileFormat == 'genbank': srcDataFileExt = 'gbff'
if options.srcFileFormat == 'embl':    srcDataFileExt = 'embl'

files = [f for f in sorted(glob.glob(gbFilesDir+'/*.'+srcDataFileExt+'.gz'))] # if os.path.isfile(f)]
for seqFile in files:

	sys.stdout.flush()
	sys.stderr.flush()
	fpIdsDNA.flush()
	fpIdsAA.flush()
	fpIdsAAadd.flush()
	fpOrgIDs.flush()

	#os.fsync(fpIdsDNA.fileno()) # ?

	if fpSeqNrDNA is not None: fpSeqNrDNA.flush()
	if fpSeqAllDNA is not None: fpSeqAllDNA.flush()
	if fpSeqAA is not None: fpSeqAA.flush()

	seqfileID = base64.b16encode(hashlib.md5(seqFile.encode('ascii')).digest())
	seqfileID=seqfileID.decode()
	print(seqfileID, seqFile),
	with gzip.open(seqFile,'rt') as fp:
		try:
			records = SeqIO.parse(fp, options.srcFileFormat) #"genbank")
		except Exception:
			continue
		if options.interpreteMode == 'A': # An assembly file, presumably containing whole genome of a single organism
			taxon = None; organism = None; assemblyID = None
			taxons = set(); organisms = set() #; Genes = dict()
			AAidPerOrgHash=set();
		for record in records:
			ok=True
			if options.doDebug: print ('\n\n\n-----Record------'); print(record)
			keywords = ''
			if 'keywords' in record.annotations:
				keywords = ';'.join(record.annotations['keywords'])
			if noFilter or record.id in SeqIdsToFilter or record.id[:-2] in SeqIdsToFilter:
				print (record.id+"\t"+record.description)
				if flagNPatch.search(keywords) is not None:
					if options.doDebug: print ('NOVEL_PATCH, locus skipped')
					continue
				if flagAltLoc.search(keywords) is not None:
					if options.doDebug: print ('ALTERNATE_LOCUS, locus skipped')
					continue
				if options.isComplete:
					if  flagWGS.search(record.description) is None \
					and flagWGS.search(keywords) is None \
					and flagComplete.search(record.description) is None \
					and flagComplete.search(keywords) is None:
						continue
				if options.interpreteMode == 'G':
					assemblyID = record.id
					taxon = None; organism = None
					taxons = set(); organisms = set() #; Genes = dict()
				gene = dict(); localgeneid = None
				CDS = dict(); CDSs = dict();
				chromosome = 'Unknown';

				#AAPTMs = list()
				source_features=OrderedDict()
				# for pasing genome assembly files
				m = re.search('(GC[FA]_\d+\.\d)', seqFile)
				if m is not None:
					assemblyID = m.group(1)
				# replace [sometimes] the one found with value from GBFF
				if record.dbxrefs is not None:
					for xrefs in record.dbxrefs:
						xref = xrefs.split(':')
						if xref[0] == 'Assembly':
							assemblyID = xref[1]

				for feature in record.features:
					if options.doDebug: print ('-----Feature------'); print(feature)
					if feature.type is not None and feature.type == 'source':
						for ft in my_source_features:
							if ft in feature.qualifiers:
								source_features[ft] = feature.qualifiers[ft][0]
						_organism = json.dumps(source_features)
						if _organism is None: continue

						if 'chromosome' in feature.qualifiers: chromosome = feature.qualifiers['chromosome'][0]
						if 'db_xref' in feature.qualifiers:
							for feature_qualifier in feature.qualifiers['db_xref']:
								if feature_qualifier[0:5] == 'taxon':
									taxons.add(int(feature_qualifier[6:])) # real tax id
						if _organism not in organisms: organisms.add(_organism)
						if options.doDebug: print ('---Source data extracted:'); print(json.dumps(feature.qualifiers))
						#if options.doDebug: print(json.dumps(source_features))

					if feature.type is not None and feature.type == 'gene':
						gene = dict()
						localgeneid = None
						if 'db_xref' in feature.qualifiers and (feature.qualifiers['db_xref'][0]).split(':')[0] == 'GeneID':
							gene['gid'] = (feature.qualifiers['db_xref'][0]).split(':')[1]
						if 'gene' in feature.qualifiers:
							gene['name'] = feature.qualifiers['gene'][0]
						if 'locus_tag' in feature.qualifiers:
							gene['locus_tag'] = feature.qualifiers['locus_tag'][0]
						if options.gidMode == 'G':
							if   'locus_tag' in gene: localgeneid = gene['locus_tag']
							elif 'gid'  in gene: localgeneid = gene['gid']
							elif 'name' in gene: localgeneid = gene['name']
						if options.gidMode == 'N':
							if   'name' in gene: localgeneid = gene['name']
							elif 'locus_tag' in gene: localgeneid = gene['locus_tag']
							elif 'gid'  in gene: localgeneid = gene['gid']
						if localgeneid is None or localgeneid == '': localgeneid = 'anonymous_gene'
						if feature.location is not None:
							gene['location'] = feature.location
							localgeneid = localgeneid + '@' + str(gene['location'])
						#Genes[localgeneid] = gene
						if options.doDebug: print ('---Gene data extracted:'); print(localgeneid); print(gene)

					if feature.type is not None and feature.type == 'CDS' and (localgeneid is not None or not options.mustHaveGeneForCDS):
						CDS = dict()
						if localgeneid is None: localgeneid='?gene?'
						if feature.location is not None:
							CDS['location'] = feature.location
						if 'product' in feature.qualifiers:
							CDS['product'] = feature.qualifiers['product'][0]
						if 'protein_id' in feature.qualifiers:
							CDS['protein_id'] = feature.qualifiers['protein_id'][0]
						if 'gene' in feature.qualifiers:
							CDS['gene'] = feature.qualifiers['gene'][0]
							if 'location' in gene:  CDS['genelocation'] = gene['location']
						if 'translation' in feature.qualifiers:
							CDS['translation'] = feature.qualifiers['translation'][0]
							CDS['proteinlength'] = len(CDS['translation'])
							if (localgeneid in CDSs):
								if (CDSs[localgeneid]['proteinlength'] < CDS['proteinlength']):
									CDSs[localgeneid] = CDS
							else: CDSs[localgeneid] = CDS
						if options.doDebug: print ('---CDS data extracted:'); print(localgeneid); print(CDS)

					#if feature.type is not None and feature.type == 'mat_peptide':
					#	print(feature.location.start,feature.location.end)
					#	print (AAseq[feature.location.start:feature.location.end+1])

				#if options.doDebug: print ('------------------all CDS-----------------'); print(CDSs)
				# end of parsing Feature table
				if options.filterTaxIdFile is not None and len(TaxIdsToFilter.intersection(taxons)) == 0 : continue
				if options.excludeTopNCBItaxon is not None and len(allExcludedTaxonsUnderTop.intersection(taxons)) > 0 : continue
				if options.topNCBItaxon is not None: taxons = allTaxonsUnderTop.intersection(taxons)
				if taxon is None:
					if len(taxons) == 0:
						taxon = None
					if len(taxons) > 0:
						taxon = str(max(taxons)) # we want the biggest NCBI tax number as it suppose to be the latest, hence the most precise tax id
					#if len(taxons) > 1:
					#	# provision for stupid cases for an assembly of a host and its pathogen
					#	sys.stderr.write(assemblyID + '\t' + ';'.join(list(map(str,taxons))) + '\n')
				if len(organisms) > 0: organism = ';'.join(organisms)


				if options.doDebug: print ('----Org IDs----'); print(seqfileOrgIDs); print(taxon); print(taxons); print(organisms)
				if organism is not None and (taxon is not None or not options.mustHaveTaxon) and (len(CDS) > 0 or not options.mustHaveCDS):
					if not noFilter:
						if   record.id      in SeqIdsToFilter: SeqIdsToFilter.remove(record.id)
						elif record.id[:-2] in SeqIdsToFilter: SeqIdsToFilter.remove(record.id[:-2])
					s = str(record.seq)
					seqlength = str(len(record.seq))
					dgstOrg = base64.b16encode(hashlib.md5((organism.replace(' ','')).encode('ascii')).digest())
					dgstOrg = dgstOrg.decode()
					dgstOrg = dgstOrg[:10]
					if options.interpreteMode == 'G': # GenBank sequence mode
						if   options.fastaFilename == 'T':
							orgID = taxon
						elif options.fastaFilename == 'D':
							orgID = dgstOrg
						elif options.fastaFilename == 'A':
							orgID = record.id
						elif options.fastaFilename == 'a':
							orgID = record.id[:-2]
						else:
							orgID = record.id
						if orgID is None: orgID = dgstOrg
						orgIDhashstr =  taxon + record.id + dgstOrg
						fastafilename = orgID
					if options.interpreteMode == 'A': # An assembly file, presumably containing whole genome of a single organism
						fastafilename = seqfileID
						orgIDhashstr  = seqfileID
						orgID         = seqfileID
						#elif options.fastaFilename.upper() == 'A' and assemblyID is not None:
						#	orgID = assemblyID
						## when assembly file contains data from several assembly versions
						#if options.fastaFilename == 'a' and assemblyID is not None:
						#	fastafilename = assemblyID.split('.')[0]

					if options.doDebug: print ('----orgIDrecord----'); print(orgID); print(taxon); print(assemblyID); print(dgstOrg); print(organism)
					orgIDrecord =  orgID + "\t" +  taxon + "\t" +  assemblyID + "\t" +  dgstOrg + "\t" +  organism + "\n"
					if orgIDhashstr not in orgIDhash:
						orgIDhash.add(orgIDhashstr)
						fpOrgIDs.write(orgIDrecord)

					# write DNA sequences
					if options.dnaFasta or options.allSeqFiles:
						dgstSeq = base64.b16encode(hashlib.md5(s.encode('ascii')).digest())
						dgstSeq = dgstSeq.decode()
						seqDNArecord =     record.id + "\t" + orgID + "\t" + chromosome + "\t" + seqlength + "\t" +  dgstSeq + "\t" + record.description + "\n"
						fastaHeader =      record.id + "\t" +  orgID
						fpIdsDNA.write(seqDNArecord)
						# per-organism files
						if options.dnaFasta:
							write_into_taxon_file(outFastaDNAdir,fastafilename,s,fastaHeader)
						# one big file all together
						if options.allSeqFiles:
							fpSeqAllDNA.write(fastaHeader + "\t" + s + "\n")
							## check id sequence is new, write to NR file if yes
							if dgstSeq not in DNAseqHash:
								DNAseqHash.add(dgstSeq)
								fpSeqNrDNA.write(dgstSeq + "\t" + s + "\n") # non-redundant DNA seqs
						#ok_fpSeqAllDNA = False

					for localgeneid,CDS in CDSs.items():
						#if options.doDebug: print(repr(CDS))
						AAid = ''; AAname = ''; AAseq = '';
						if 'translation' in CDS:
							AAseq = CDS['translation']
							#seqlength = str(len(AAseq))
							if 'protein_id' in CDS:
								AAid = CDS['protein_id']
							if 'product' in CDS:
								AAname = CDS['product']
							if AAid > '':
								if options.excludeAAIDdups: # filter out duplicating protein ids within this organism
									if AAid in AAidPerOrgHash:
										continue
									else: AAidPerOrgHash.add(AAid)
								if options.fetchFromENTREZ and 'polyprotein' in AAname.lower():
									while True:
										try:
											cnx = Entrez.efetch(db="protein", id=AAid, rettype="gb", retmode="xml")
											record_pp = Entrez.read(cnx)
											cnx.close()
											break
										except Exception:
											sys.stderr.write ("\nRestarting... "+ AAid +"\n")

									if options.doDebug: print ('-----Entrez------'); print(record_pp)
									#sys.stderr.write(json.dumps(record_pp[0]['GBSeq_feature-table'], sort_keys=True, indent=4))
									GBFeature_keys=set()
									pprotein_interval_annotation=dict()
									if 'GBSeq_sequence' in record_pp[0]:
										sequence = record_pp[0]['GBSeq_sequence'].upper()
									else: continue
									for rrr in record_pp[0]['GBSeq_feature-table']: # find all keys
										GBFeature_keys.add(rrr['GBFeature_key'])
									if 'mat_peptide' in GBFeature_keys: the_pp_feature = 'mat_peptide'
									else: the_pp_feature = 'Region'
									for rrr in record_pp[0]['GBSeq_feature-table']:
										if options.doDebug: print ('-----Entrez feature ------'); print(rrr)
										if rrr.get('GBFeature_key') is not None and rrr['GBFeature_key'] == the_pp_feature:
											#sys.stderr.write(json.dumps(rrr, sort_keys=True, indent=4))
											AAregSeq = None; AAregID = None
											pprotein_interval_annotation.clear()
											#pprotein_interval_annotation.update(protein_annotation)
											pprotein_interval_annotation['region_from_to'] = rrr['GBFeature_location']
											pprotein_interval_annotation.update(get_gb_gualifiers(rrr['GBFeature_quals'],polyprotein_cv[the_pp_feature]))
											if the_pp_feature ==  'mat_peptide' and 'peptide' in pprotein_interval_annotation:
												AAregSeq = pprotein_interval_annotation.pop('peptide')
											elif the_pp_feature ==  'Region' or the_pp_feature ==  'mat_peptide':
												if 'product' in pprotein_interval_annotation:
													AAname = pprotein_interval_annotation['product']
												elif 'region_name' in pprotein_interval_annotation:
													AAname = pprotein_interval_annotation['region_name']
												elif 'name' in pprotein_interval_annotation:
													AAname = pprotein_interval_annotation['name']
												#if 'note' in pprotein_interval_annotation and pprotein_interval_annotation.get('note') is not None:
												#	pprotein_interval_annotation['product'] += '; '+pprotein_interval_annotation.pop('note')
												seq_first = int(rrr['GBFeature_intervals'][0]['GBInterval_from'])
												seq_last  = int(rrr['GBFeature_intervals'][0]['GBInterval_to'])
												AAregSeq = sequence[seq_first-1:seq_last]
												if AAregSeq is not None:
													dgstSeq = base64.b16encode(hashlib.md5(AAregSeq.encode('ascii')).digest())
													dgstSeq = dgstSeq.decode()
													if 'protein_id' in pprotein_interval_annotation:
														AAregID = pprotein_interval_annotation.get('protein_id')
													else:
														AAregID = AAname+"@"+str(seq_first)+":"+str(seq_last)+"@"+AAid
											if AAregSeq is not None and AAregID is not None :
												fastaHeader = AAregID  + "\t" +  orgID
												seqAArecord = AAregID + "\t" + record.id + "\t" + orgID + "\t" +  localgeneid + "@" +  chromosome   + "\t" + dgstSeq + "\t" + AAname
												fpIdsAA.write(seqAArecord + "\n")
												if options.aaFasta:
													write_into_taxon_file(outFastaAAdir,fastafilename,AAregSeq,fastaHeader)
								else:
									if options.aaFasta:
										dgstSeq = base64.b16encode(hashlib.md5(AAseq.encode('ascii')).digest())
										dgstSeq = dgstSeq.decode()
										fastaHeader =      AAid + "\t" + orgID
										seqAArecord =      AAid + "\t" + record.id + "\t" + orgID + "\t" +  localgeneid + "@" +  chromosome   + "\t" + dgstSeq + "\t" + AAname
										fpIdsAA.write(seqAArecord + "\n")
										write_into_taxon_file(outFastaAAdir,fastafilename,AAseq,fastaHeader)

						elif 'protein_id' in CDS: # those protein ids that are only mentioned in CDS record without providing their AA seq
							fpIdsAAadd.write(CDS['protein_id']+ "\n")
				else:
					if options.doDebug: print ('--- NO DATA WRITTEN ---')
#cursor.close
if not noFilter and len(SeqIdsToFilter) > 0:
	for id in SeqIdsToFilter:
		fpAbsentIDs.write(id + "\n")
