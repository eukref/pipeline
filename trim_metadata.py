import os
import sys
import glob
import argparse

parser = argparse.ArgumentParser(
    description='Trim metadata table down to the taxa in clustered reference tree')
parser.add_argument(
    '-i',
    '--input_fasta',
    help='Clustered Reference Fasta File',
    required=True)
parser.add_argument(
    '-m',
    '--input_metadatafile',
    help='Metadata File',
    required=True)
parser.add_argument(
    '-o',
    '--output_trimmed_metadatafile',
    help='Trimmed Metadata File Output',
    required=True)
    
args = parser.parse_args()

fas_infile = args.input_fasta
meta_file = args.input_metadatafile
out_file = args.output_trimmed_metadatafile

infile = open(fas_infile)
line = infile.read()
infile.close()
seqs = line[1:].split('\n>')

centroids_accs = []

for seq in seqs:
	acc = seq.split()[0]
	acc = acc.split('_')[0]
	if acc != 'OUTGROUP':
		centroids_accs.append(acc)
 

infile = open(meta_file)
lines = infile.readlines()
infile.close()

out = open(out_file, "w")
out.write(lines[0])

counter = 0

for line in lines[1:]:
	if line.split()[0] in centroids_accs:
		out.write(line)
		counter = counter + 1

if counter == len(centroids_accs):
	print 'Metadata file trimmed to %s records' % (len(centroids_accs))
else:
	print 'WARNING: There is %s records recovered from Metadata file and %s sequences in fasta file, check that no records are missing' % (len(centroids_accs), counter) 
 		
out.close()

