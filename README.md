# epinorm-genbank
Parser for GenBank files

typical example to parse gzipped GBFF files from subdir "src" without creating fasta sequences:
parseGBFF.py --sourcefeatures=OSITHGDC --interpretemode=G src > 1.log 2>1.err.log &

A run creates several files in the current directory. Important are standard output [redirected to 1.log] and id.org.tab containing data with NCBI organism tax id in the 2nd column and organism description in JSON-like format in the last column.

Plesae, use --help option for more options
