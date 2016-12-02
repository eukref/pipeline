#!/usr/bin/env python

print "\nScript for parsing GenBank metadata and renaming fasta file for annotation." 
#print "Contributors: Javier del Campo and Laura Wegener Parfrey"
#print "22 November 2016"
print "\nrun: python eukref_gbmetadata.py -h for help.\n"

import os
import argparse
import re
from Bio import SeqIO
from Bio.Blast import NCBIXML




parser = argparse.ArgumentParser(
    description='Rename fasta sequences and create metadata text file for FigTree annotation')
parser.add_argument(
    '-gb',
    '--input_gb_file_fp',
    help='Path to GenBank formatted file with accessions to be renamed. Get this file from the NCBI nucleotide database online (see pipeline overview)',
    required=False)
parser.add_argument(
    '-t', 
    '--ref_tax',
    help='Path to reference database taxonomy file formatted accession \t taxonomy. E.g. SILVA 128 full taxa map file. Reference database taxonomy will be added to the metadata file',
    required=False)
parser.add_argument(
    '-i',
    '--input_fasta_file_fp',
    help='Path to Fasta file to be renamed. May be current_DB.fas, current_DB.clustered.fas, or other. Header MUST either be 1) in standard GenBank format (e.g. >gi|ginumber|gb|accession| ) or 2) begin with the accession number. Get this file from the NCBI nucleotide database online (see pipeline overview)',
    required=False)
parser.add_argument(
    '--outgroup',
    help='Path to Fasta file with outgroups.',
    required=False)
parser.add_argument(
    '-o',
    '--output_fasta_file_fp',
    help='Path to output fasta file',
    required=False)
parser.add_argument(
    '-m',
    '--output_metadata_fp',
    help='Output metadata file in tab delimited format',
    required=False)

args = parser.parse_args()

ref_tax = args.ref_tax
input_gb_file_fp = args.input_gb_file_fp
input_fasta_file_fp = args.input_fasta_file_fp
outgroup = args.outgroup
output_fasta_file_fp = args.output_fasta_file_fp
output_metadata_fp = args.output_metadata_fp


outfasta = open(output_fasta_file_fp, "w")
outmeta = open(output_metadata_fp, "w")

##############################
#### function definitions ####
##############################

# function just loops through outgroup fasta file and adds OUTGROUP to line. writes to outfile. 
def rename_outgroup(infile, outfile):
	for line in open(infile, "U"):
		if line.startswith(">"):
			line = line.replace(">" , ">OUTGROUP__")
			outfile.write(line)
		else:
			outfile.write(line)



# function to print metadata in tab delimited format based on gb formatted input file. 
# version if SILVA or PR2 reference taxonomy file passed
def metadata_retrieve_ref(infile, outfile, ref_accessions):
	# original script from here
	OUT = outmeta
	OUT.write("accession\tgenbank_taxonomy\treference_taxonomy\torganism_name_gb\tclone_name_gb\tsource_gb\tenvironment_gb\thost_gb\tcountry_gb\tpublication_gb\tauthors_gb\tjournal_gb\n")
	result_handle = open(infile, "U")
	# array to make sure each accession is uniq. 
	uniq_acc = []
	gbfiles = SeqIO.parse(result_handle, 'gb')
	for rec in gbfiles:
		acc = rec.id
		# strip off .1 from accessions 
		clean_acc = re.sub(r'\.[1-9]', '', acc)
		#if already have seen accession move onto next record. Else append accession to uniq_acc list and 
		if clean_acc in uniq_acc:
			next
		else:
			uniq_acc.append(clean_acc)
			#default = 'NA'
			#reference_taxonomy = ref_accessions.get('clean_acc', default)
			if clean_acc in ref_accessions:
				reference_taxonomy = ref_accessions[clean_acc]
			else:
				reference_taxonomy = 'NA'
			source = rec.features[0]
			if 'taxonomy' in rec.annotations:	
				taxonomy = ";".join(rec.annotations['taxonomy'])
			if 'organism' in rec.annotations:
				organism = rec.annotations['organism']
			else:
				organism = "NA"
			if 'clone' in source.qualifiers:
				clone = source.qualifiers['clone'][0]
			else:
				clone = "NA"
			if 'environmental_sample' in source.qualifiers:
				environmental_sample = "Environmental"
			else:
				environmental_sample = "Isolate"
			if 'isolation_source' in source.qualifiers:
				isolation_source = source.qualifiers['isolation_source'][0]
			else:
				isolation_source = "NA"
			if 'host' in source.qualifiers:
				host = source.qualifiers['host'][0]
			else:
				host = "NA"
			if 'country' in source.qualifiers:
				country = source.qualifiers['country'][0]
			else:
				country = "NA"
			if 'references' in rec.annotations:
				pubref = rec.annotations['references'][0]
				authors = pubref.authors
				title = pubref.title
				journal = pubref.journal
			else:
				title = "NA"
				authors = "NA"
				journal = "NA"
			fields = [clean_acc, taxonomy, reference_taxonomy, organism, clone, environmental_sample, isolation_source, host, country, title, authors, journal]
			OUT.write("\t".join(fields)+ "\n")
	OUT.close()		

# function to print metadata in tab delimited format based on gb formatted input file. 
def metadata_retrieve(infile, outfile):
	accessions = {}
	# original script from here
	OUT = outmeta
	OUT.write("accession\tgenbank_taxonomy\torganism_name_gb\tclone_name_gb\tsource_gb\tenvironment_gb\thost_gb\tcountry_gb\tpublication_gb\tauthors_gb\tjournal_gb\n")
	result_handle = open(infile, "U")
	uniq_acc = []
	gbfiles = SeqIO.parse(result_handle, 'gb')
	for rec in gbfiles:
		acc = rec.id
		# need to make sure each accession is unique. 
		# strip off .1 from accessions 
		clean_acc = re.sub(r'\.[1-9]', '', acc)
		#if already have seen accession move onto next record. Else append accession to uniq_acc list and 
		if clean_acc in uniq_acc:
			next
		else:
			uniq_acc.append(clean_acc)
			source = rec.features[0]
			if 'taxonomy' in rec.annotations:	
				taxonomy = ";".join(rec.annotations['taxonomy'])
			if 'organism' in rec.annotations:
				organism = rec.annotations['organism']
			else:
				organism = "NA"
			if 'clone' in source.qualifiers:
				clone = source.qualifiers['clone'][0]
			else:
				clone = "NA"
			if 'environmental_sample' in source.qualifiers:
				environmental_sample = "Environmental"
			else:
				environmental_sample = "Isolate"
			if 'isolation_source' in source.qualifiers:
				isolation_source = source.qualifiers['isolation_source'][0]
			else:
				isolation_source = "NA"
			if 'host' in source.qualifiers:
				host = source.qualifiers['host'][0]
			else:
				host = "NA"
			if 'country' in source.qualifiers:
				country = source.qualifiers['country'][0]
			else:
				country = "NA"
			if 'references' in rec.annotations:
				pubref = rec.annotations['references'][0]
				authors = pubref.authors
				title = pubref.title
				journal = pubref.journal
			else:
				title = "NA"
				authors = "NA"
				journal = "NA"
			fields = [clean_acc, taxonomy, organism, clone, environmental_sample, isolation_source, host, country, title, authors, journal]
			OUT.write("\t".join(fields)+ "\n")
	OUT.close()			

def rename_sequences_ref(in_gb, infile, ref_acc, outfile):
	# initialize gb_tax dict to store accession, last taxonomy level, and organism name
	gb_tax = {}
	#result_handle = open(in_gb, "U")
	# parse gb file; block to go through gb file for acc not in ref
	gbfiles = SeqIO.parse(in_gb, 'gb')
	# make dict of acc, last level, organism
	for rec in gbfiles:
		acc = rec.id
		clean_acc = re.sub(r'\.[1-9]', '', acc)
		#print clean_acc
		if clean_acc in ref_acc:
			next
		else:
			if 'organism' in rec.annotations:
				organism = rec.annotations['organism']
			if 'taxonomy' in rec.annotations:	
				# get last level of taxonomy
				taxonomy = rec.annotations['taxonomy']
				# add accession plus last taxonomy level and organism name to tax dict
				name = taxonomy[-1]+"_"+ organism
				name = name.replace(" ", "_")
				#print name
				gb_tax[clean_acc] = name

	# go through fasta file. Match accessions first to reference taxonomy (ref_acc) then to gb_tax
	for line in open(infile, "U"):
		if line.startswith(">"):
			#seq_acc = None
			if "_" in line:
				seq_acc = line.split("_")[0]
				seq_acc = seq_acc.replace(">", "")
			# case of old gb format with >gi|noginumber|gb|KT210044 
			elif "|" in line:
				acc = seq.split('|')[3]
				seq_acc = re.sub(r'\.[1-9]', '', acc)
				seq_acc = seq_acc.replace(">", "")
			# case of just accession # in header, or accession followed by white space
			else:
				seq_acc = line.split()[0]
				seq_acc = seq_acc.replace(">", "")
			if seq_acc in ref_acc:
				taxonomy = ref_acc[seq_acc]
				#name = None
				if "; __" in taxonomy:
					last_level = taxonomy.split("; __")
					name = '_'.join(last_level[-2:])
					#tax[seq_acc] = '_'.join(last_level[-2:-1])
				else:
					last_level = taxonomy.split(";")
					name = '_'.join(last_level[-2:])
				outfile.write(">%s_%s\n" % (seq_acc,name))
			elif seq_acc in gb_tax:
				outfile.write(">%s_%s\n" % (seq_acc, gb_tax[seq_acc]))
			# incase accession is missing. 
			else:
				outfile.write("%s\n" % (line.strip()))
		else:
			outfile.write("%s\n" % (line.strip()))
				


# ###################################################################
# ######################### SCRIPT ITSELF ###########################
# ###################################################################

# Make dict of silva taxonomy (or PR2) if provided. 
# how to handle SILVA? read whole taxonomy file into dictionary with key, value as accession, taxonomy? 
# Maybe here only want to make a dict of accession and whole taxonomy. do not split unless necessary. 
ref_accessions = {}
if args.ref_tax is not None:
	for line in open(ref_tax, "U"):
		# split by \t
		acc = line.strip().split('\t')
		ref_accessions[acc[0]] = acc[1]
		#print ref_accessions[acc[0]]			


# run  eukref_gbmetadata.py version that also reports reference taxonomy if given gb file AND given reference taxonomy file (e.g. Silva taxonomy file)
if args.input_gb_file_fp is not None and args.ref_tax is not None:
	# accessions is a dictionary of accessions in file (with name?)
	metadata_retrieve_ref(input_gb_file_fp, outmeta, ref_accessions)

#run eukref_gbmetadata.py if gb file given. 
if args.input_gb_file_fp is not None and args.ref_tax is None:
	metadata_retrieve(input_gb_file_fp, outmeta)
	print "metadata file is %s" % (outmeta)
	

if args.outgroup is not None:
	# run outgroup renaming and write outgroup seqs to outfile
	rename_outgroup(outgroup, outfasta)
	

# doesn't really make sense to only use ref info - should also use gb record. 
if args.input_fasta_file_fp is not None and args.ref_tax is not None and args.input_gb_file_fp is not None: 
	# need to make this function
	rename_sequences_ref(input_gb_file_fp, input_fasta_file_fp, ref_accessions, outfasta)
	print "fasta file for tree is %s" % (outfasta)
	
outfasta.close()
outmeta.close()
