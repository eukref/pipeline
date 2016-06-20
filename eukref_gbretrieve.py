import os
import pickle
import argparse

parser = argparse.ArgumentParser(
    description='Iterative Taxon Level Specific Blast')
parser.add_argument(
    '-i',
    '--db',
    help='Dataset of sarting sequences in fasta format, inputy filename',
    required=True)
parser.add_argument(
    '-dbnt',
    '--path_to_nt',
    help='Path to NCBI nt database, string',
    required=True)
parser.add_argument(
    '-dbsi',
    '--path_to_silva',
    help='Path to Silva database, formated as DNA, string',
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
    '--percent_identity_silva',
    help='percent identity for blast hits against SILVA, float',
    required=True)
args = parser.parse_args()

db = args.db
path_to_nt = args.path_to_nt
path_to_silva = args.path_to_silva
num_seq_per_round = args.num_seqs_per_round
cpu = args.cpu
group_name = args.group_name
percent_id_nt = float(args.percent_identity_nt)
percent_id_silva = float(args.percent_identity_silva)
blast_method = args.nt_blast_method

print db
print path_to_nt
print path_to_silva
print num_seq_per_round
print cpu
print group_name

# loads in the tax dictionary
tax_d = pickle.load(open('tax_d.bin', 'rb'))

Renaming_d = {}


# runs blast against NT database
def nt_blast(query_file, num_seqs, outf):
    os.system(
        "blastn -task %s -query %s -db %s -num_threads %s -max_target_seqs %s -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames scomnames sblastnames sskingdoms' -out %s" %
        (blast_method,
         query_file,
         path_to_nt,
         str(cpu),
         str(num_seqs),
         outf))
    infile = open(outf)
    lines = infile.readlines()
    infile.close()
    acc = []
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
    #os.system('blastn -task megablast -query %s -db %s -num_threads %s -max_target_seqs %s -outfmt 6 -out %s' % (query_file, path_to_silva, cpu, 1, outf))
    os.system(
        './usearch -usearch_local %s -db %s -strand both -maxhits 1 -id 0.8 -threads %s -blast6out %s' %
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
            print jaba
            if tax_d[ac_hit].count('Bacteria') > 0 or tax_d[
                    ac_hit].count('Archaea') > 0:
                bad_hits[line.split('|')[0]] = 'This is bad hit'
                print ac_hit, 'bacteria'
            # elif tax_d[ac_hit].count('uncultured') and tax_d[ac_hit].count(group_name) == 0:
            #	print 'confusing'

            #	confusing_sequences.append(line.split('|')[3])
            elif float(p_id) > 95.0 and tax_d[ac_hit].count(group_name) == 0 and tax_d[ac_hit] not in allowed:
                print 'bad hit'
                bad_hits[line.split('|')[3]] = 'This is bad hit'
                print tax_d[ac_hit]
                print ac_hit
            elif float(p_id) < float(percent_id_silva):
                bad_hits[line.split('|')[3]] = 'This is bad hit'
                print 'KURVA?'
                print ac_hit
                print float(p_id)
            else:
                coord_dict[line.split('|')[3]] = (
                    line.split()[6], line.split()[7])
                acc.append(line.split('|')[3])
    print acc
    print len(acc)
    return (acc, coord_dict)
    print confusing_sequences
    print len(bad_hits)


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


# removes sequences shorter then length
def trim_short_seqs(fnamein, fnameout, length):
    d = load_fasta_file(fnamein)
    out = open(fnameout, 'w')
    print d
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
    print d
    return d


def run_uchime(fnamein, fnameout):
    print fnamein
    whynot = open(fnamein, 'r')
    line = whynot.read()
    print line
    whynot.close()
    d = {}
    seqs = line[1:].split('\n>')
    print seqs
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
    if os.path.isfile('usearch') == False:
        os.system('usearch -uchime_ref temp_chime_in.fas -db gb203_pr2_all_10_28_97p_noorg.udb -nonchimeras temp_no_chimeras.fas -chimeras temp_chimeras.fas -strand plus')
    else:
        os.system('./usearch -uchime_ref temp_chime_in.fas -db gb203_pr2_all_10_28_97p_noorg.udb -nonchimeras temp_no_chimeras.fas -chimeras temp_chimeras.fas -strand plus')

    infile = open('temp_no_chimeras.fas')
    line = infile.read()
    infile.close()
    seqs = line[1:].split('\n>')
    out = open(fnameout, 'w')
    print seqs
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


def rename_sequences(infile, outfile):
    infile = open(infile)
    line = infile.read()
    infile.close()
    seqs = line[1:].split('\n>')
    out = open(outfile, 'w')
    for seq in seqs:
        try:
            new_name = seq.split('|')[3] + '|' + '_'.join(Renaming_d[seq.split('|')[3]][0].split()) + '|' + '_'.join(
                Renaming_d[seq.split('|')[3]][1].split()) + '|' + '_'.join(Renaming_d[seq.split('|')[3]][2].split())
        except KeyError:
            new_name = seq.split('\n')[0]
        out.write('>%s\n%s\n' % (new_name, ''.join(seq.split('\n')[1:])))
    out.close()


#####################################################################
############################SCRIPT ITSELF############################
#####################################################################
#os.system('makeblastdb -in %s -out %s -dbtype nucl' % (path_to_silva, path_to_silva))
shortout = open('short_seqs_acs.txt', 'w')
chimeraout = open('chimera_seqs_acs.txt', 'w')
bad_hits = {}
allowed = []
confusing_sequences = []

#########################Analyze initial DB##########################

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
        out.write('>gi|noginumber|gb|noaccesssion%s|%s\n' % (acc_count, seq))
        acc_count = acc_count + 1
out.close()

os.system('cp %s %s_old_version' % (db, db))
os.system('cp %s_reformated %s' % (db, db))
os.system('cp %s current_DB.fas' % (db))

print 'Runing uchime on initial fasta file'
run_uchime(db, 'temp_initial_db.fas')
print 'Removing short sequences from initial fasta file'
trim_short_seqs('temp_initial_db.fas', 'new_round.fas', 500)

# counter of cycles
c = 0

db_dict = load_fasta_file(db)

#########################Start Cycling DB##########################
# It is 'yes' until
keep_running = 'yes'

while keep_running == 'yes':
    c = c + 1
    print 'starting cycle %s' % (int(c))
#	run usearch
    print 'running usearch'
    if os.path.isfile('usearch') == False:
        os.system(
            'usearch -sortbylength new_round.fas -fastaout DB.sorted.fas -minseqlength 64')
    else:
        os.system(
            './usearch -sortbylength new_round.fas -fastaout DB.sorted.fas -minseqlength 64')

    if os.path.isfile('usearch') == False:
        os.system(
            'usearch -cluster_smallmem DB.sorted.fas -id 0.97 -centroids temp_DB_clust.fas -uc temp_DB_clust.uc')
    else:
        os.system(
            './usearch -cluster_smallmem DB.sorted.fas -id 0.97 -centroids temp_DB_clust.fas -uc temp_DB_clust.uc')

#	remove seqs shorter 500
#	trim_short_seqs('temp_DB_clust.fas', 'temp_DB_clust_500.fas', 500)
#	print 'removing short sequences'
#	remove species that are already in the database
    print 'running blast'
    hits_acs = nt_blast(
        'temp_DB_clust.fas',
        num_seq_per_round,
        'temp_DB_clust_500.blastntab')
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

# Kill cycle if no new sequences
    if len(hits_acs) == 0:
        print 'NO NEW SEQUENCES'
        keep_running = 'no'
# Keep running next cycle
    else:
        #		write accessions into text file
        out = open('temp_to_get_acs.txt', 'w')
        for i in hits_acs:
            out.write('%s\n' % (i))
        out.close()
#		recover fasta file for accessions
        os.system(
            'blastdbcmd -entry_batch temp_to_get_acs.txt -db %s -out temp_results.fas' %
            (path_to_nt))
        trim_short_seqs('temp_results.fas', 'temp_results_500.fas', 500)
# blast against Silva remove <70%ID and >95%ID other groups, This resturn
# accession numbers
        print 'blasting against SILVA'
        silva_blast_results = silva_blast(
            'temp_results_500.fas', 'temp_results.silvatab')
        silva_parsed_acs = silva_blast_results[0]
        silva_coord_dict = silva_blast_results[1]
        print silva_coord_dict
#		Write the accession numbers into text file
        if silva_parsed_acs != []:
            out = open('temp_to_get_acs.txt', 'w')
            for i in silva_parsed_acs:
                out.write('%s\n' % (i))
            out.close()
#			pull passed sequences from the NCBI db
            os.system(
                'blastdbcmd -entry_batch temp_to_get_acs.txt -db %s -out temp_results_silva_untrimmed.fas' %
                (path_to_nt))
#			run uchime
#			Trim sequences based on best hit
            infile = open('temp_results_silva_untrimmed.fas')
            line = infile.read()
            infile.close()
            out = open('temp_results_silva.fas', 'w')
            seqs = line[1:].split('\n>')
            for i in seqs:
                new_seq = ''.join(i.split('\n')[1:])
                hit_len = int(silva_coord_dict[i.split('|')[3]][
                              0]) - int(silva_coord_dict[i.split('|')[3]][1])
                print hit_len
                if hit_len > 0:

                    if float(hit_len) / float(len(new_seq)) < 0.8:
                        new_seq = new_seq[int(silva_coord_dict[i.split('|')[3]][0]):int(
                            silva_coord_dict[i.split('|')[3]][1])]
                    else:
                        new_seq = new_seq
                elif hit_len < 0:
                    print -float(hit_len) / float(len(new_seq))
                    if -float(hit_len) / float(len(new_seq)) < 0.8:
                        print new_seq
                        new_seq = new_seq[int(silva_coord_dict[i.split('|')[3]][0]):int(
                            silva_coord_dict[i.split('|')[3]][1])]
                        print new_seq
                    else:
                        new_seq = new_seq
                out.write('>%s\n%s\n' % (i.split('\n')[0], new_seq))
            out.close()
            print 'removing chimeras'
            print 'SHITSHITSHIT'
            prdel = open('temp_results_silva.fas')
            print prdel
            print prdel.read()
            run_uchime(
                'temp_results_silva.fas',
                'temp_results_silva_uchime.fas')
            prdel = open('temp_results_silva.fas')
            print prdel.readlines()
#			remove short sequences
            trim_short_seqs(
                'temp_results_silva_uchime.fas',
                'temp_results_silva_uchime_500.fas',
                500)
#			merge new sequences with old database
            os.system(
                'cat current_DB.fas temp_results_silva_uchime_500.fas > new_DB.fas')
            os.system('mv new_DB.fas current_DB.fas')
#			restart cycle with new sequences
            os.system('mv temp_results_silva_uchime_500.fas new_round.fas')
#			remove temporary files - move to ne folder
            os.system('mkdir cycle_%s' % (c))
            os.system('mv temp* cycle_%s' % (c))
            os.system('cp current_DB.fas cycle_%s' % (c))
#			check there are new sequences left
            infile = open('new_round.fas')
            line = infile.read()
            infile.close()
            if line == '':
                keep_running = 'no'
        else:
            keep_running = 'no'


print confusing_sequences
rename_sequences('current_DB.fas', 'current_DB_done.fas')

#########TO ADD##########
# re-recover short seqs
# re-recover chimera seqs
# deal with Genome sequences
# add gi restricted blast to remove the already removed seqs - not it is
# little clumsy
