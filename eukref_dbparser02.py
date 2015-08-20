'''
usage
python annotateclusters.py thecamoebids.clustered.fasttree.txt thecamoebids.clusters.uc expandedclusters.txt

'''


import os
import sys
import glob
import re

infile = open(sys.argv[1],"r")  ### tab delimited outfile of annotatetree.py
lines = infile.readlines()
infile.close()
clusterfile = open(sys.argv[2], "r") ### *.uc file
clusterlines = clusterfile.readlines()
clusterfile.close()
outfile = open(sys.argv[3],"w") ### tab delimited outfile


for line in clusterlines:
	cluster = line.split("\t")[9]
	cluster = cluster.strip()
	fullname = line.split("\t")[8]
	
	if fullname.startswith("gi+"):
		accession = line.split("+")[3]

	if fullname.startswith("gi|"):
		accession = line.split("|")[3]

	else:
		pass

	if cluster == "*":
		for i in lines:
			if i.startswith(accession):
				#print i
				outfile.write(i)
			else:
				pass	

	else:
		for i in lines:
			if i.startswith(cluster):
				annotation = i.split("\t")[1]
				outfile.write(accession + "\t" + annotation)
			else:
				pass

outfile.close()
sys.exit()
