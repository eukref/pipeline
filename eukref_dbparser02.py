'''
usage
python annotateclusters.py thecamoebids.clustered.fasttree.txt thecamoebids.clusters.uc expandedclusters.txt

'''

import sys

infile = open(sys.argv[1],"r")  ### tab delimited outfile of annotatetree.py
lines = infile.readlines()
infile.close()
clusterfile = open(sys.argv[2], "r") ### *.uc file
clusterlines = clusterfile.readlines()
clusterfile.close()
outfile = open('temp_expandedcluster.txt',"w") ### tab delimited outfile

for line in clusterlines:
	#line = line.replace("+","|")
	cluster = line.split("\t")[9]
	cluster = cluster.strip()
	fullname = line.split("\t")[8]
	if fullname.count("noginumber") > 0:
		accession = line.split('|')[-1]
		accession = accession.split()[0]

	elif fullname.startswith("gi+"):
		accession = fullname.split('+')[3]

	elif fullname.startswith("gi|"):
		accession = fullname.split("|")[3]

	
	else:
		pass
        accession.strip()
#        print accession
	#print fullname
	#print accession
	if cluster == "*":
		for i in lines:
#                        print i
			if i.startswith(accession):
				#print i
				outfile.write(i)
			else:
				pass	

	else:
		for i in lines:
			if cluster.count(i.split()[0]) > 0:
				annotation = i.split("\t")[1]
				outfile.write(accession + "\t" + annotation)
			else:
				pass

outfile.close()

#remove duplicated lines
infile = open('temp_expandedcluster.txt')
lines = infile.readlines()
infile.close()

list_lines = []

for line in lines:
  list_lines.append(line.strip())

new_lines = set(list_lines)
out = open(sys.argv[3],'w')

for i in new_lines:
	out.write('%s\n' % (i))

out.close()

sys.exit()
