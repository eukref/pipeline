#!/usr/bin/env python
print "\nScript to pull FigTree annotation into tab delimited metadata file." 
print "Contributors: Laura Wegener Parfrey and Evan Morien"
print "November 2016"
print "\nrun: python eukref_treeparse.py -h for help.\n"


import argparse

parser = argparse.ArgumentParser(
    description='Pull FigTree annotation into tab delimited metadata file')
parser.add_argument(
    '-t',
    '--annotated_tree_file_path',
    help='Path to annotated tree from FigTree. Add annotations to "taxa" in FigTree and do file "Save as" to save annotations to .tre file',
    required=True)
parser.add_argument(
    '-m',
    '--input_metadata_file',
    help='Path to tab delimited metadata file.',
    required=True)

args = parser.parse_args()

tree = args.annotated_tree_file_path
metadata = args.input_metadata_file

outfile = metadata.strip(".txt")+"_out.txt"
print "Metadata with annotation added to last column is in file %s\n" % (outfile)
out_metadata = open(outfile, "w")
#print out_metadata
#open(out_metadata, "w")

tree_annotations = {}
for line in open(tree, "U"):
	if "[&!" in line and not ";" in line:
		# get annotation
		annotation = line.split('"')[1]
		accession = line.split('[&!')[0]
		accession.strip()
		accession = accession.split('_')[0]
		accession = accession.split('.')[0]
		accession = accession.split()[0]
		accession = accession.strip("'")
		tree_annotations[accession] = annotation
		#print "%s\t%s" % (accession, annotation)
		
for line in open(metadata, "U"):
		
	parts = line.strip().split("\t")
	acc = parts[0]
	if line.startswith("Accession"):
		parts.append("Annotation")
	elif acc in tree_annotations:
		parts.append(tree_annotations[acc])
	else:
		parts.append("none")
	new_line = "\t".join(parts)
	out_metadata.write("%s\n" % (new_line))
		
out_metadata.close()
