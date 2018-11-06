#!/usr/bin/env python2

"""
Writen by Serafim Nenarokov.
6.11.2018

run: ./eukref_recover_headers.py
"""
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def extract_acc(header):
    acc = header.split()[0]
    return acc.split('.')[0]


def main():
    cur_db_path = 'current_DB.fas'
    cur_db_final_path = 'current_DB_final.fas'
    acc_path = 'accessions_current_DB.txt'

    with open(cur_db_path) as cur_db_f, open(cur_db_final_path, 'w') as cur_db_final_f:
        with open(acc_path, 'w') as acc_f:

            records = {}

            for rec in SeqIO.parse(cur_db_f, 'fasta'):

                if '>' in rec.description:
                    headers = rec.description.split(' >')
                else:
                    headers = [rec.description]

                for header in headers:
                    acc = extract_acc(header)

                    new_rec = SeqRecord(rec.seq, id=acc,
                                        name=acc,
                                        description=header)

                    if acc not in records:
                        records[acc] = new_rec

            for acc, rec in records.iteritems():
                cur_db_final_f.write(rec.format('fasta'))
                acc_f.write(extract_acc(acc) + "\n")

    # check the number of lines
    cur_db_cnt = int(subprocess.check_output("grep -c '>' %s" % cur_db_path, shell=True).strip())
    cur_db_final_cnt = int(subprocess.check_output("grep -c '>' %s" % cur_db_final_path, shell=True).strip())
    acc_cnt = int(subprocess.check_output("cat %s | wc -l" % acc_path, shell=True).strip())

    if cur_db_cnt == cur_db_final_cnt == acc_cnt:
        print("Everything is fixed correctly (%s seqs)" % acc_cnt)
    else:
        err_msg = "ERROR: (current_DB: %s records, current_DB_final: %s records, accessions: %s lines)"
        print(err_msg % (cur_db_cnt, cur_db_final_cnt, acc_cnt))


if __name__ == '__main__':
    main()
