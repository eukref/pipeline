#!/usr/bin/env python
print "\nexpands clusters" 
print "Contributors: Matthew Brown, Laura Wegener Parfrey, and Evan Morien"
print "November 2016"
print "\nrun: python eukref_dbparser02.py -h for help.\n"
print "\nExample run: python eukref_dbparser02.py -i Metadata_File -r reference_metadata_file -c current_DB.clusters.uc  -o expandedclusters.txt\n"


import sys
import argparse

parser = argparse.ArgumentParser(
    description='expand clusters')
parser.add_argument(
    '-i',
    '--input_metadata',
    help='metadata',
    required=True)
parser.add_argument(
    '-r',
    '--input_reference_metadata',
    help='reference_metadata',
    required=True)
parser.add_argument(
    '-c',
    '--input_clusterfile',
    help='cluster file ==> *.uc file',
    required=True)
parser.add_argument(
    '-o',
    '--output_expandedclusters',
    help='metadata file expanded to all sequences in reference database with annotations',
    required=True)
    
args = parser.parse_args()

meta_infile = args.input_metadata
trim_file = args.input_reference_metadata
out_file = args.output_expandedclusters
cluster_file = args.input_clusterfile
    
infile = open(cluster_file)
lines = infile.readlines()
infile.close()

centroid_d = {}

for line in lines:
	if line.split()[0] == 'S':
		accession_rec = line.split('\t')[8]
		print accession_rec
		accession_rec = accession_rec.split()[0]
                print accession_rec
		accession_rec = accession_rec.split("_")[0]
                print accession_rec
		accession_rec = accession_rec.split(".")[0]
		print accession_rec
		centroid_d[line.split('\t')[1]] = accession_rec
	else:
		pass
print centroid_d
acc_centroid_d = {}
for line in lines:
#	print line
#	print centroid_d[line.split('\t')[1]]
#	print acc_centroid_d
	accession_rec = line.split('\t')[8]
        accession_rec = accession_rec.split()[0]
      	accession_rec = accession_rec.split("_")[0]
       	accession_rec = accession_rec.split(".")[0]
	acc_centroid_d[accession_rec] = centroid_d[line.split('\t')[1]]

#print acc_centroid_d	
#print acc_centroid_d['AF315604']

infile = open(meta_infile)
lines = infile.readlines()
infile.close()
number_col = len(lines[0].split('\t'))
meta_d = {}
for i in lines[1:]:
	meta_d[i.split('\t')[0]] = i.strip()
	
infile = open(trim_file)
lines = infile.readlines()
infile.close()

trim_d = {}

for line in lines[1:]:
#	print line
#	print number_col
	annotation = '\t'.join( line.split('\t')[number_col:])
	trim_d[line.split('\t')[0]] = annotation.strip()
print trim_d
print meta_d
print acc_centroid_d
out = open(out_file, 'w')

for acc in meta_d:
#	print acc
	try:
		out.write('%s\t%s\n' % (meta_d[acc], trim_d[acc_centroid_d[acc]]))
	except KeyError:
		print 'WARNING accession number %s does no exist in cluster file' % (acc)
		out.write('%s\t%s\n' % (meta_d[acc], trim_d[acc]))
out.close()


print acc_centroid_d
	
 









	




