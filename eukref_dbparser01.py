"""

usage
python annotatetree.py thecamoebids.clustered.fasttree.nex thecamoebids.clustered.fasttree.txt

"""

import sys
import re

infile = open(sys.argv[1], "r")  # annotated tree nexus file
lines = infile.readlines()
infile.close()
outfile = open(sys.argv[2], "w")  # tab delimited outfile

for line in lines:
    if re.search('!name=', line):
        print line
        annotation = line.split('!name="')[1]
        annotation = annotation.split('"')[0]
        if not line.startswith("\t'gi+") and not line.startswith("\t'gi|"):
            pass
        else:
            if line.startswith("\t'gi+"):
                separator = "+"
            else:
                separator = "|"            
            accession = line.split(separator)[3]
            accession = accession.split(separator)[0]
            if accession == "noaccesssion":
                accession2 = line.split(separator)[4]
                accession2 = accession2.split("'")[0]
                outfile.write(accession2 + "\t" + annotation + "\n")
            else:
                outfile.write(accession + "\t" + annotation + "\n")
outfile.close()

sys.exit()
