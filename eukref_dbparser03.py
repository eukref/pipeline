"""
usage
removes sequences that are annotated as "remove" in your figtree and passed through annotatetree.py then annotateclusters.py
python editfasta.py expandedclusters.txt current_DB.fas current_DB-cleared.fas
"""

import sys


def oneline(infasta):
    infile = open('%s' % (infasta), "r")
    lines = infile.readlines()
    infile.close()
    outfile = open('%s' % (infasta), 'w')
    for i, line in enumerate(lines):
        if line[0] == ('>'):
            if i > 0:
                outfile.write("\n")
            outfile.write(line)
        else:
            line = line.strip()
            outfile.write(line)
    outfile.write("\n")
    outfile.close()

# tab delimited outfile of annotateclusters.py ==> "expandedclusters.txt"
infile = open(sys.argv[1], "r")
lines = infile.readlines()
infile.close()

infasta = sys.argv[2]  # big fasta file ==> current_DB.fas
oneline(infasta)
infastafile = open(infasta, "r")
fastalines = infastafile.readlines()
infastafile.close()

# cleared fastafile of "remove" sequences ==> "current_DB-cleared.fas"
outfile = open(sys.argv[3], "w")

removallist = []

for line in lines:
    removal = line.split("\t")[1]
    removal = removal.strip()
    id = line.split("\t")[0]
    if removal == ("remove" or "REMOVE" or "Remove"):
        removallist.append(id)
print removallist
for j in range(len(fastalines)):
    line = fastalines[j]
    if line[0] == ">":

        if line.count("noginumber") > 0:
            accession = line.split("|")[-1]
            accession = accession.strip()
            print accession
            if accession not in removallist:
                print 'Not_deleteing'
                outfile.write(line)
                next_l = fastalines[j + 1]
                outfile.write(next_l)

        elif line.startswith(">gi+"):
            accession = line.split("+")[3]
            if accession not in removallist:
                outfile.write(line)
                next_l = fastalines[j + 1]
                outfile.write(next_l)

        elif line.startswith(">gi|"):
            accession = line.split("|")[3]
            if accession not in removallist:
                outfile.write(line)
                next_l = fastalines[j + 1]
                outfile.write(next_l)

        else:
            print 'sssssss'
            outfile.write(line)
            next_l = fastalines[j + 1]
            outfile.write(next_l)

    else:
        pass

outfile.close()
sys.exit()
