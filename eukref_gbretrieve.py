"""
Writen by Martin Kolisko.
29.8.2016

run: python eukref_gbretrieve.py -h for usage.

This pipeline requires USEARCH (32bit version is OK) and NCBI Blast to be installed. If you are using linux and cannot
install programs as administrator, you can install both under your account and add them to your PATH or put executables
into the folder you are using this script from. 

You will also need localy placed NCBI NT database, which can be downloaded from NCBI ftp server.
"""



import os
import pickle
import argparse
import re

print "Genbank SSU rRNA gene parser." 
print "Contributors: Martin Kolisko, Javier del Campo, Laura Wegener Parfrey and Frederick Mahe"
print "29 August 2016"
print "\n\nrun: python eukref_gbretrieve.py -h for help.\n\nThis pipeline requires USEARCH (32bit version is OK) and NCBI Blast to be installed. If you are using linux and cannot\ninstall programs as administrator, you can install both under your account and add them to your PATH.  \n\nYou will also need localy placed NCBI NT database, which can be downloaded from NCBI ftp server."
print "\n\n"
parser = argparse.ArgumentParser(
    description='Iterative Taxon Level Specific Blast')
parser.add_argument(
    '-i',
    '--db',
    help='Dataset of sarting sequences in fasta format. MUST be in standard GenBank format (>gi|ginumber|gb|accession| ), or have the accession number followed by a space (>accession ), input filename',
    required=True)
parser.add_argument(
    '-dbnt',
    '--path_to_nt',
    help='Path to NCBI nt database, string',
    required=True)
parser.add_argument(
    '-dbsi',
    '--path_to_ReferenceDB',
    help='Path to Reference database, formated as DNA, string',
    required=True)
parser.add_argument(
    '-n',
    '--num_seqs_per_round',
    help='Number of sequences recovered per blast, integer',
    required=True)
parser.add_argument(
    '-p',
    '--cpu',
    help='Number of cpu cores used for blast, integer',
    required=True)
parser.add_argument(
    '-g',
    '--group_name',
    help='Name of taxon group to be sampled',
    required=True)
parser.add_argument(
    '-m',
    '--nt_blast_method',
    help='"megablast" or "blastn"',
    required=True)
parser.add_argument(
    '-idnt',
    '--percent_identity_nt',
    help='percent identity for blast hits against NT, float',
    required=True)
parser.add_argument(
    '-idsi',
    '--percent_identity_ReferenceDB',
    help='percent identity for blast hits against Reference DB, float',
    required=True)
parser.add_argument(
    '-td',
    '--taxonomy_dict',
    help='path to taxonomy dictionary - file tax_d.bin',
    required=True)
    
args = parser.parse_args()

db = args.db
path_to_nt = args.path_to_nt
path_to_silva = args.path_to_ReferenceDB
num_seq_per_round = args.num_seqs_per_round
cpu = args.cpu
group_name = args.group_name
percent_id_nt = float(args.percent_identity_nt)
percent_id_silva = float(args.percent_identity_ReferenceDB)
blast_method = args.nt_blast_method
taxonomy_dict_file = args.taxonomy_dict


print db
print path_to_nt
print path_to_silva
print num_seq_per_round
print cpu
print group_name

# loads in the tax dictionary
tax_d = pickle.load(open('tax_d.bin', 'rb'))

Renaming_d = {}

os.system("usearch -makeudb_usearch Reference_DB.fas -output Reference_DB.udb")
os.system("usearch -makeudb_usearch gb203_pr2_all_10_28_97p_noorg.fasta -output gb203_pr2_all_10_28_97p_noorg.udb")

# runs blast against NT database
def nt_blast(query_file, num_seqs, outf):
    OUTFORMAT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sblastnames sskingdoms"
    os.system(
        "blastn -task %s -query %s -db %s -num_threads %s -max_target_seqs %s -outfmt '%s' -out %s" %
        (blast_method,
         query_file,
         path_to_nt,
         str(cpu),
         str(num_seqs),
         OUTFORMAT,
         outf))
    infile = open(outf)
    lines = infile.readlines()
    infile.close()
    acc = []
    # this is the part that renames the seq.  
    for line in lines:
        if float(line.split()[2]) >= float(percent_id_nt):
            hit_column = line.split()[1]
            acc.append(hit_column.split('|')[3])
            Renaming_d[hit_column.split('|')[3]] = [line.split(
                '\t')[-4], line.split('\t')[-2], line.split('\t')[-1]]
    accset = set(acc)

    return accset


def silva_blast(query_file, outf):
    coord_dict = {}
    os.system(
        'usearch -usearch_local %s -db %s -strand both -maxhits 1 -id 0.8 -threads %s -blast6out %s' %
        (query_file, path_to_silva, cpu, outf))
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
            if tax_d[ac_hit].count('Bacteria') > 0 or tax_d[
                    ac_hit].count('Archaea') > 0:
                bad_hits[line.split('|')[0]] = 'This is bad hit'
             #   print ac_hit, 'bacteria'
            elif float(p_id) > 95.0 and tax_d[ac_hit].count(group_name) == 0 and tax_d[ac_hit] not in allowed:
             #   print 'bad hit'
                bad_hits[line.split('|')[3]] = 'This is bad hit'
             #   print tax_d[ac_hit]
             #   print ac_hit
            elif float(p_id) < float(percent_id_silva):
                bad_hits[line.split('|')[3]] = 'This is bad hit'
             #   print ac_hit
              #  print float(p_id)
            else:
                coord_dict[line.split('|')[3]] = (
                    line.split()[6], line.split()[7])
                acc.append(line.split('|')[3])
#    print acc
#    print len(acc)
    return (acc, coord_dict)
 #   print confusing_sequences
 #   print len(bad_hits)


# load fasta file into dictionary
def load_fasta_file(fname):
    infile = open(fname)
    line = infile.read()
    infile.close()
    seqs = line[1:].split('\n>')
    d = {}
    for seq in seqs:
        d[seq.split('\n')[0]] = ''.join(seq.split('\n')[1:])
    return d


# removes sequences shorter than 500 bp in length. At end of run prints accessions of short seqs to short_seqs_acc.txt
def trim_short_seqs(fnamein, fnameout, length):
    d = load_fasta_file(fnamein)
    out = open(fnameout, 'w')
    #print d
    for i in d:
        if i != '':
            if len(d[i]) > int(length) and len(d[i]) < 5000:
                out.write('>%s\n%s\n' % (i, d[i]))
            else:
                shortout.write('%s\n' % (i.split('|')[3]))


# makes_empt_dict of accession numbers in fasta file
def make_ac_dict(fname):
    d = {}
    infile = open(fname)
    line = infile.read()
    infile.close()
    seqs = line[1:].split('\n>')
    for seq in seqs:
        d[seq.split('|')[3]] = ''
    #print d
    return d


# run uchime to search for chimeric seqs. 
def run_uchime(fnamein, fnameout):
    whynot = open(fnamein, 'r')
    line = whynot.read()
    whynot.close()
    d = {}
    seqs = line[1:].split('\n>')
    #print seqs
    out = open('temp_chime_in.fas', 'w')
    for seq in seqs:
        d[seq.split('|')[3]] = seq
        out.write(
            '>%s;size=1;\n%s\n' %
            (seq.split('|')[3],
             ''.join(
                seq.split('\n')[
                    1:])))
    out.close()
    os.system('usearch -uchime_ref temp_chime_in.fas -db %s -nonchimeras temp_no_chimeras.fas -chimeras temp_chimeras.fas -strand plus' % (path_to_silva))

    infile = open('temp_no_chimeras.fas')
    line = infile.read()
    infile.close()
    seqs = line[1:].split('\n>')
    out = open(fnameout, 'w')
    #print seqs
    if line.count('>') != 0:
        for seq in seqs:
            out.write('>%s\n' % (d[seq.split(';')[0]]))
    else:
        pass

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
    out.close()


# function for outputting a list of all accessions to be used further. 
def out_accessions(infile, out_fasta_file, out_accessions_file):
	out_acc = open(out_accessions_file, 'w')
	out_fasta = open(out_fasta_file, 'w')
	accessions = []
	for line in open(infile, "U"):
		if line.startswith(">"):
			line = line.split(">")[1]
			try:
				# find the accession from GenBank format
				acc = line.split('|')[3]
				clean_acc = re.sub(r'\.[1-9]', '', acc)
				if clean_acc in accessions:
					next
				else:
					accessions.append(clean_acc)
					out_acc.write("%s\n" % (clean_acc))
					out_fasta.write(">%s\n" % (clean_acc))
			except IndexError:
			# if not in old genbank format assume in form of >accession.1 other info, with a space delimiter
				acc = seq.split()[0]
				clean_acc = re.sub(r'\.[1-9]', '', acc)
				if clean_acc in accessions:
					next
				else:
					accessions.append(clean_acc)				
					out_acc.write("%s\n" % (clean_acc))
					out_fasta.write(">%s\n" % (clean_acc))
		else:
			out_fasta.write("%s\n" % (line.strip()))
				
			
#function for renaming sequences to place in tree. Start with header default from GenBank. End up with >accession_name_from_genbank
# def rename_sequences(infile, outfile):
# 	out = open(outfile, 'w')
# 	for l in open(infile, "U"):
# 		if l.startswith(">"):
# 			line = l.split(">")[1]  
# 			# bunch of lines that remove annoying bits from fasta header
#  			line = line.strip().replace(" 18S ribosomal RNA gene, partial sequence", "")
#  			line = line.replace(" 18S ribosomal RNA, partial sequence", "")
#  			line = line.replace(" 18S small subunit ribosomal RNA_gene, complete sequence", "")
#  			line = line.replace(" gene for 18S rRNA, partial sequence", "")
#  			line = line.replace(" partial 18S rRNA gene", "")
#  			line = line.replace(" small subunit ribosomal RNA gene, partial sequence", "")
#  			acc = line.split('|')[3].strip('.[1-9]')
#  			label = line.split(" ")
#  			new_label = "_".join(label[1:])
#  			out.write(">%s_%s\n" % (acc, new_label))
#  		else:
#  			out.write("%\n" % (l.strip()))
#  	out.close()

# Martin's function for renaming sequences
# def rename_sequences(infile, outfile):
#     infile = open(infile)
#     line = infile.read()
#     infile.close()
#     seqs = line[1:].split('\n>')
#     out = open(outfile, 'w')
#     for seq in seqs:
#         try:
#             new_name = seq.split('|')[3] + '|' + '_'.join(Renaming_d[seq.split('|')[3]][0].split()) + '|' + '_'.join(
#                 Renaming_d[seq.split('|')[3]][1].split()) + '|' + '_'.join(Renaming_d[seq.split('|')[3]][2].split())
#         except KeyError:
#             new_name = seq.split('\n')[0]
#         out.write('>%s\n%s\n' % (new_name, ''.join(seq.split('\n')[1:])))
#     out.close()


# ###################################################################
# ######################### SCRIPT ITSELF ###########################
# ###################################################################
shortout = open('short_seqs_acc.txt', 'w')
chimeraout = open('chimera_seqs_acc.txt', 'w')
bad_hits = {}
allowed = []
confusing_sequences = []


# ###################### Analyze initial DB #########################
# Format New Sequences
infile = open(db)
line = infile.read()
infile.close()
acc_count = 0
out = open('%s_reformated' % (db), 'w')
seqs = line[1:].split('\n>')
for seq in seqs:
    try:
        ac = seq.split('|')[3]
        out.write('>%s\n' % (seq))
    except IndexError:
    	# if not in old genbank format assume in form of >accession.1 other info, with a space delimiter
		accession = seq.split(" ")[0]
		out.write('>gi|noginumber|gb|%s| %s\n' % (accession, seq))
    acc_count = acc_count + 1
out.close()

os.system('cp %s %s_old_version' % (db, db))
#os.system('cp %s_reformated %s' % (db, db))
os.system('cp %s_reformated current_DB.fas' % (db))

print 'Runing uchime on initial fasta file'
run_uchime('current_DB.fas', 'temp_initial_db.fas')
print 'Removing short sequences from initial fasta file'
trim_short_seqs('temp_initial_db.fas', 'new_round.fas', 500)

# counter of cycles
c = 0

db_dict = load_fasta_file(db)

# ###################### Start Cycling DB #########################
while True:
    c = c + 1
    print 'starting cycle %s' % (int(c))
    # run usearch
    print 'running usearch'
    os.system('usearch -sortbylength new_round.fas -fastaout DB.sorted.fas -minseqlength 64')
    os.system('usearch -cluster_smallmem DB.sorted.fas -id 0.97 -centroids temp_DB_clust.fas -uc temp_DB_clust.uc')

    # remove seqs shorter 500
    # trim_short_seqs('temp_DB_clust.fas', 'temp_DB_clust_500.fas', 500)
    # print 'removing short sequences'
    # remove species that are already in the database
    print 'running blast'
    hits_acs = nt_blast(
        'temp_DB_clust.fas',
        num_seq_per_round,
        'temp_DB_clust_500.blastntab')
    # remove_repeats when compared to current_DB.fas
    dict_existing = make_ac_dict('current_DB.fas')
    hits_acs2 = list(hits_acs)
    for i in hits_acs2:
        try:
            yup = dict_existing[i]
            hits_acs.remove(i)
        except KeyError:
            pass
    #print len(hits_acs)

    # remove repeated BAD hits
    for i in hits_acs2:
        try:
            yup = bad_hits[i]
            hits_acs.remove(i)
        except KeyError:
            pass
    #print len(hits_acs)

    # Kill cycle if no new sequences
    if len(hits_acs) == 0:
        print 'NO NEW SEQUENCES'
        break
    # Keep running next cycle
    else:
        # write accessions into text file
        out = open('temp_to_get_acs.txt', 'w')
        for i in hits_acs:
            out.write('%s\n' % (i))
        out.close()
        # recover fasta file for accessions
        os.system(
            'blastdbcmd -entry_batch temp_to_get_acs.txt -db %s -out temp_results.fas' %
            (path_to_nt))
        trim_short_seqs('temp_results.fas', 'temp_results_500.fas', 500)
        # blast against Silva remove <70%ID and >95%ID other groups, This returns
        # accession numbers
        print 'blasting against SILVA'
        silva_blast_results = silva_blast(
            'temp_results_500.fas', 'temp_results.silvatab')
        silva_parsed_acs = silva_blast_results[0]
        silva_coord_dict = silva_blast_results[1]
        #print silva_coord_dict
        # Write the accession numbers into text file
        if silva_parsed_acs != []:
            out = open('temp_to_get_acs.txt', 'w')
            for i in silva_parsed_acs:
                out.write('%s\n' % (i))
            out.close()
            # pull passed sequences from the NCBI db
            os.system(
                'blastdbcmd -entry_batch temp_to_get_acs.txt -db %s -out temp_results_silva_untrimmed.fas' %
                (path_to_nt))
            # run uchime
            # Trim sequences based on best hit
            infile = open('temp_results_silva_untrimmed.fas')
            line = infile.read()
            infile.close()
            out = open('temp_results_silva.fas', 'w')
            seqs = line[1:].split('\n>')
            for i in seqs:
                new_seq = ''.join(i.split('\n')[1:])
                hit_len = int(silva_coord_dict[i.split('|')[3]][
                              0]) - int(silva_coord_dict[i.split('|')[3]][1])
                #print hit_len
                if hit_len > 0:

                    if float(hit_len) / float(len(new_seq)) < 0.8:
                        new_seq = new_seq[int(silva_coord_dict[i.split('|')[3]][0]):int(
                            silva_coord_dict[i.split('|')[3]][1])]
                    else:
                        new_seq = new_seq
                elif hit_len < 0:
                    #print -float(hit_len) / float(len(new_seq))
                    if -float(hit_len) / float(len(new_seq)) < 0.8:
                        #print new_seq
                        new_seq = new_seq[int(silva_coord_dict[i.split('|')[3]][0]):int(
                            silva_coord_dict[i.split('|')[3]][1])]
                        #print new_seq
                    else:
                        new_seq = new_seq
                out.write('>%s\n%s\n' % (i.split('\n')[0], new_seq))
            out.close()

            prdel = open('temp_results_silva.fas')
            run_uchime(
                'temp_results_silva.fas',
                'temp_results_silva_uchime.fas')
            prdel = open('temp_results_silva.fas')
            # remove short sequences
            trim_short_seqs(
                'temp_results_silva_uchime.fas',
                'temp_results_silva_uchime_500.fas',
                500)
            # merge new sequences with old database
            os.system(
                'cat current_DB.fas temp_results_silva_uchime_500.fas > new_DB.fas')
            os.system('mv new_DB.fas current_DB.fas')
            # restart cycle with new sequences
            os.system('mv temp_results_silva_uchime_500.fas new_round.fas')
            # remove temporary files - move to ne folder
            os.system('mkdir cycle_%s' % (c))
            os.system('mv temp* cycle_%s' % (c))
            os.system('cp current_DB.fas cycle_%s' % (c))
            # check there are new sequences left
            infile = open('new_round.fas')
            line = infile.read()
            infile.close()
            if line == '':
                break
        else:
            break

# run function to rename sequences in output fasta file
#rename_sequences('current_DB.fas', 'current_DB_tree_names.fas')

# run function to output a list of accessions in the final DB
out_accessions('current_DB.fas', 'current_DB_final.fas', 'accessions_current_DB.txt')

# ####### TO ADD #########
# re-recover short seqs
# re-recover chimera seqs
# deal with Genome sequences
# add gi restricted blast to remove the already removed seqs - not it is
# little clumsy
