#!/usr/bin/env python

import os
import pickle
import argparse
from sys import argv
import re
import argparse
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
    '--ref_tax',
    help='Path to reference database taxonomy file formatted accession \t taxonomy. E.g. SIVLA 128 full taxa map file. Reference database taxonomy will be added to the metadata file',
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


# this function should 1) Check input fasta for duplicates. 2) write out accessions and include name? or taxonomy? or SILVA taxonomy?

# function for renaming fasta file header using metadata information. 
#Should add in processing of the gb file?
# out sequence (which should be viewable in tree) name should be >accession_taxonomy_name
def rename_sequences(infile, outfile):
	for line in open(infile, 'U'):
		if line.startswith(">"):
			seq_header = line.split(">")[1]
			seq_header = seq_header.strip().replace(" 18S ribosomal RNA gene, partial sequence", "")
			seq_header = seq_header.replace(" 18S ribosomal RNA, partial sequence", "")
			seq_header = seq_header.replace(" 18S small subunit ribosomal RNA_gene, complete sequence", "")
			seq_header = seq_header.replace(" gene for 18S rRNA, partial sequence", "")
			seq_header = seq_header.replace(" partial 18S rRNA gene", "")
			seq_header = seq_header.replace(" small subunit ribosomal RNA gene, partial sequence", "")
			acc = line.split('|')[3]
			clean_acc = re.sub(r'\.[1-9]', '', acc)
			label = seq_header.split(" ")
			new_label = "_".join(label[1:])
			outfile.write(">%s_%s\n" % (clean_acc, new_label.strip()))
		else:
			outfile.write("%s\n" % (line.strip()))


# function to print metadata in tab delimited format based on gb formatted input file. 
# version if SILVA or PR2 reference taxonomy file passed
def metadata_retrieve_ref(infile, outfile, ref_accessions):
	# original script from here
	OUT = outmeta
	OUT.write("Accession\tTaxonomy\tReference_taxonomy\tOrganism\tclone\tSource\tEnvironment\tHost\tCountry\tPublication\tAuthors\tJournal\n")
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
	OUT.write("Accession\tTaxonomy\tOrganism\tclone\tSource\tEnvironment\tHost\tCountry\tPublication\tAuthors\tJournal\n")
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
	print "retrieving metadata"


if args.outgroup is not None:
	# run outgroup renaming and write outgroup seqs to outfile
	rename_outgroup(outgroup, outfasta)
	

# if passed a fasta file as input will run function to rename fasta headers
if args.input_fasta_file_fp is not None:
	rename_sequences(input_fasta_file_fp, outfasta)

outfasta.close()
outmeta.close()
