#!/usr/bin/env python
# Author: Javier del Campo
# Date: August 2016
# Usage: python eukref_gbmetadata.py outfile.txt infile.gb
# Requires BioPython.  Download and installation instructions: http://biopython.org/wiki/Download

import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
OUT = open(sys.argv[1], 'w')
OUT.write("Accession\tTaxonomy\tName\tSource\tEnvironment\tHost\tCountry\tPublication\tAuthors\tJournal\n")
result_handle = open(sys.argv[2])
gbfiles = SeqIO.parse(result_handle, 'gb')
for rec in gbfiles:
    acc = rec.id
    source = rec.features[0]
    if 'taxonomy' in rec.annotations:	
        taxonomy = ";".join(rec.annotations['taxonomy'])
    if 'clone' in source.qualifiers:
        clone = source.qualifiers['clone'][0]
    elif 'organism' in rec.annotations:
        clone = rec.annotations['organism']
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
    fields = [acc, taxonomy, clone, environmental_sample, isolation_source, host, country, title, authors, journal]
    OUT.write("\t".join(fields)+ "\n")
OUT.close()	
