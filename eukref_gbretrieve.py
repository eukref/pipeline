import os
import sys
import glob
import pickle
import argparse

parser = argparse.ArgumentParser(description='Iterative Taxon Level Specific Blast')
parser.add_argument('-i', '--db', help='Dataset of sarting sequences in fasta format, inputy filename', required=True)
parser.add_argument('-dbnt', '--path_to_nt', help='Path to NCBI nt database, string', required=True)
parser.add_argument('-dbsi', '--path_to_silva', help='Path to Silva database, formated as DNA, string', required=True)
parser.add_argument('-n', '--num_seqs_per_round', help='Number of sequences recovered per blast, integer', required=True)
parser.add_argument('-p', '--cpu', help='Number of cpu cores used for blast, integer', required=True)
parser.add_argument('-g', '--group_name', help='Name of taxon group to be sampled', required=True)
parser.add_argument('-m', '--nt_blast_method', help='"megablast" or "blastn"', required=True)
parser.add_argument('-id', '--percent_identity', help='percent identity for blast hits, float', required=True)
args = parser.parse_args()

db = args.db
path_to_nt = args.path_to_nt
path_to_silva = args.path_to_silva
num_seq_per_round = args.num_seqs_per_round
cpu = args.cpu
group_name = args.group_name
percent_id = args.percent_identity
blast_method = args.nt_blast_method

print db
print path_to_nt
print path_to_silva
print num_seq_per_round
print cpu
print group_name

shortout = open('short_seqs_acs.txt', 'w')
chimeraout = open('chimera_seqs_acs.txt', 'w')
bad_hits = {}

#loads in the tax dictionary
tax_d = pickle.load(open('ac_taxonomy.dat','rb'))

#runs blast against NT database
def nt_blast(query_file, num_seqs, outf):
	os.system('blastn -task %s -query %s -db %s -num_threads %s -max_target_seqs %s -outfmt 6 -out %s' % (blast_method, query_file, path_to_nt, str(cpu), str(num_seqs), outf))
	infile = open(outf)
	lines = infile.readlines()
	infile.close()
	acc = []
	for line in lines:
		if float(line.split()[2]) >= float(percent_id):
			hit_column = line.split()[1]
			acc.append(hit_column.split('|')[3])
	accset = set(acc)

	return accset

def silva_blast(query_file, outf):
	os.system('blastn -task megablast -query %s -db %s -num_threads %s -max_target_seqs %s -outfmt 6 -out %s' % (query_file, path_to_silva, cpu, 1, outf))
	infile = open(outf)
	lines = infile.readlines()
	infile.close()
	done = []
	acc = []
	for line in lines:
		p_id = line.split()[2]
		ac_hit = line.split()[1]
		ac_hit = ac_hit.split('.')[0]
		jaba = 0 
		if ac_hit not in done:
			done.append(line.split('|')[3])
			jaba = jaba + 1
			print jaba
			if float(p_id) > 95.0 and tax_d[ac_hit].count(group_name) == 0:
				bad_hits[line.split('|')[3]] = 'This is bad hit'
				print tax_d[ac_hit]
				print ac_hit
			elif float(p_id) < float(percent_id):
				bad_hits[line.split('|')[3]] = 'This is bad hit'
				print 'KURVA?'
				print ac_hit
				print float(p_id)
			else:
				acc.append(line.split('|')[3])
	print acc
	print len(acc)
	return acc
	

#load fasta file into dictionary
def load_fasta_file(fname):
	infile = open(fname)
	line = infile.read()
	infile.close()
	seqs = line[1:].split('\n>')
	d = {}
	for seq in seqs:
		d[seq.split('\n')[0]] = ''.join(seq.split('\n')[1:])
	return d

#removes sequences shorter then length
def trim_short_seqs(fnamein, fnameout, length):
	d = load_fasta_file(fnamein)
	out = open(fnameout,'w')
	for i in d:
		if len(d[i]) > int(length) and len(d[i]) < 5000:
			out.write('>%s\n%s\n' % (i, d[i]))
		else:
			shortout.write('%s\n' % (i.split('|')[3]))

#makes_empt_dict of accession numbers in fasta file
def make_ac_dict(fname):
	d = {}
	infile = open(fname)
	line = infile.read()
	infile.close()
	seqs = line[1:].split('\n>')
	for seq in seqs:
		d[seq.split('|')[3]] = ''
	print d
	return d

def run_uchime(fnamein, fnameout):
	infile = open(fnamein)
	line = infile.read()
	infile.close()
	d = {}
	seqs = line[1:].split('\n>')
	out = open('temp_chime_in.fas','w')
	for seq in seqs:
		d[seq.split('|')[3]] = seq
		out.write('>%s;size=1;\n%s\n' % (seq.split('|')[3], ''.join(seq.split('\n')[1:])))
	out.close()
	os.system('./usearch8.0.1623_i86linux32 -uchime_denovo temp_chime_in.fas -nonchimeras temp_no_chimeras.fas -chimeras temp_chimeras.fas')
	infile = open('temp_no_chimeras.fas')
	line = infile.read()
	infile.close()
	seqs = line[1:].split('\n>')
	out = open(fnameout,'w')
	for seq in seqs:
		out.write('>%s\n' %  (d[seq.split(';')[0]]))
	
	try:
		infile = open('temp_chimera.fas')
		line = infile.read()
		infile.close()
		seqs = line[1:].split('\n>')
		for seq in seqs:
			chimeraout.write('%s\n' % (seq.split('|')[3]))
		os.system('rm temp_chime_in.fas temp_chimeras.fas temp_no_chimeras.fas')
	except IOError:
		pass
	
	

	

#####################################################################
############################SCRIPT ITSELF############################	
#####################################################################


#########################Analyze initial DB##########################
print 'Runing uchime on initial fasta file'
run_uchime(db, 'temp_initial_db.fas')
print 'Removing short sequences from initial fasta file'
trim_short_seqs('temp_initial_db.fas', 'new_round.fas', 500)


#counter of cycles
c = 0


db_dict = load_fasta_file(db)

#It is 'yes' until 
keep_running = 'yes'

while keep_running == 'yes':
	c = c + 1
	print 'starting cycle %s' % (int(c))
#	run usearch
	print 'running usearch'
	os.system('./usearch -sortbylength new_round.fas -fastaout DB.sorted.fas -minseqlength 64')
	os.system('./usearch -cluster_smallmem DB.sorted.fas -id 0.97 -centroids temp_DB_clust.fas -uc temp_DB_clust.uc')
#	remove seqs shorter 500
#	trim_short_seqs('temp_DB_clust.fas', 'temp_DB_clust_500.fas', 500)
#	print 'removing short sequences'
#	remove species that are already in the database
	print 'running blast'
	hits_acs = nt_blast('temp_DB_clust.fas', num_seq_per_round, 'temp_DB_clust_500.blastntab')
#	remove_repeats when compared to current_DB.fas
	dict_existing = make_ac_dict('current_DB.fas')
	hits_acs2 = list(hits_acs)
	for i in hits_acs2:
		try:
			yup = dict_existing[i]
			hits_acs.remove(i)
		except KeyError:
			pass
	print len(hits_acs)


#	remove repeated BAD hits	
	for i in hits_acs2:
		try:
			yup = bad_hits[i]
			hits_acs.remove(i)
		except KeyError:
			pass
	print len(hits_acs)
	
	

#Kill cycle if no new sequences
	if len(hits_acs) == 0:
		print 'NO NEW SEQUENCES'
		keep_running = 'no'
#Keep running next cycle
	else:
#		write accessions into tect file
		out = open('temp_to_get_acs.txt','w')
		for i in hits_acs:
			out.write('%s\n' % (i))
		out.close()
#		recover fasta file for accessions
		os.system('blastdbcmd -entry_batch temp_to_get_acs.txt -db %s -out temp_results.fas' % (path_to_nt))
#		blast against Silva remove <70%ID and >95%ID other groups, This resturn accession numbers
		print 'blasting against SILVA'
		silva_parsed_acs = silva_blast('temp_results.fas', 'temp_results.silvatab')
#		Write the accession numbers into text file
		out = open('temp_to_get_acs.txt','w')
		for i in silva_parsed_acs:
			out.write('%s\n' % (i))
		out.close()
#		pull passed sequences from the NCBI db
		os.system('blastdbcmd -entry_batch temp_to_get_acs.txt -db %s -out temp_results_silva.fas' % (path_to_nt))
#		run uchime
		print 'removing chimeras'
		run_uchime('temp_results_silva.fas', 'temp_results_silva_uchime.fas')
#		remove short sequences
		trim_short_seqs('temp_results_silva_uchime.fas', 'temp_results_silva_uchime_500.fas', 500)
#		merge new sequences with old database
		os.system('cat current_DB.fas temp_results_silva_uchime_500.fas > new_DB.fas')
		os.system('mv new_DB.fas current_DB.fas')
#		restart cycle with new sequences
		os.system('mv temp_results_silva_uchime_500.fas new_round.fas')
#		remove temporary files - move to ne folder
		os.system('mkdir cycle_%s' % (c))
		os.system('mv temp* cycle_%s' % (c))
#		check there are new sequences left
		infile = open('new_round.fas')
		line = infile.read()
		infile.close()
		if line == '':
			keep_running = 'no'



#########TO ADD##########
#re-recover short seqs
#re-recover chimera seqs
#deal with Genome sequences
#add gi restricted blast to remove the already removed seqs
