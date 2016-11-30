#!/usr/bin/env python
print "\nScript to clean fasta headers to include only accession numbers. Input in fasta file wiht headers in old genbank format." 
print "Contributors: Laura Wegener Parfrey"
print "November 2016"
print "\nrun: python clean_fasta.py -h for help.\n"


import argparse

parser = argparse.ArgumentParser(
    description='Clean fasta headers')
parser.add_argument(
    '-i',
    '--input_fasta_file_fp',
    help='Path to input fasta file. Header must have accession number first followed by white space or a _. E.g. ">Accession  " or ">Accession_". May be current_DB.fas, current_DB.clustered.fas, or current_DB_final.fas or other. ',
    required=True)
parser.add_argument(
    '-s', 
    '--list_of_accessions',
    help='Path to file with list of accessions, or accessions in first column of tab delimited file.',
    required=True)
parser.add_argument(
    '-o',
    '--output_fasta_file_fp',
    help='Path to output fasta file',
    required=True)


args = parser.parse_args()

list_of_accessions = args.list_of_accessions
input_fasta_file_fp = args.input_fasta_file_fp
output_fasta_file_fp = args.output_fasta_file_fp

# make list of accessions
acc_remove = []
for line in open(list_of_accessions, "U"):
	if "_" in line:
		acc = line.split('_')[0]
		acc_remove.append(acc)
		next
	else:
		# in case no _. will split on any whitespace
		acc = line.split()[0]
		acc_remove.append(acc)

#process fasta file
# open input fasta file and read whole thing into infile
infile = open(input_fasta_file_fp).read()
# get seqs
seqs = infile.split('\n>')
out = open(output_fasta_file_fp, 'w')
for seq in seqs:
	# get accession
	header = seq.split()
	accession = None
	# isolate accession. With underscore or alone. 
	if "_" in header[0]:
		accession = header[0].split('_')[0]
		accession = accession.strip('>')
	else:
		accession = header[0].strip('>')
	# skip if accession in remove list
	if accession in acc_remove:
		next
	else:
		out.write('>%s\n' % (seq))
print "fasta file with sequences removed is %s" % (output_fasta_file_fp)		
out.close()