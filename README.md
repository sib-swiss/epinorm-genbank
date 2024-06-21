# epinorm-genbank
Parser for GenBank files, requires BioPython.
Here is a typical example to parse gzipped GBFF files from subdir "src" without creating fasta sequences:

parseGBFF.py --sourcefeatures=OSITHGDC --interpretemode=G src > 1.log 2>1.err.log &

A run creates several files in the current directory. Important are the standard output [redirected to 1.log], and "id.org.tab" containing data with NCBI organism tax id in the 2nd column and organism description in JSON-like format in the last column.
Please, use --help option for more options


After the above step, the resulting file id.org.tab might be additionally normalized using parseParseGBFF.py, using a command like:

parseParseGBFF.py -1 id.org.tab NCBItaxNames.csv > TBE.GB.csv

Please see --help for description
