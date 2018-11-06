#!/usr/bin/env python2

"""
Writen by Serafim Nenarokov.
6.11.2018

run: ./eukref_recover_headers.py -i <path to original fasta>
                                 -o <path to fixed fasta>
                                 -acc_list <path to accession list>
"""

import argparse

from Bio import SeqIO


def prepare_args():
    description = 'Recover the original fasta headers'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input",
                        required=True,
                        type=argparse.FileType('r'),
                        help="Path to the input current_db.fasta file.")

    parser.add_argument("-o", "--output",
                        required=True,
                        type=argparse.FileType('w'),
                        help="Path to the output fasta file.")

    parser.add_argument("--acc_list",
                        required=True,
                        type=argparse.FileType('w'),
                        help="Path where accession list will be saved.")

    return parser.parse_args()


def extract_acc(header):
    return header.split()[0]


def main():
    opts = prepare_args()

    with opts.input as input_f, opts.output as out_f, opts.acc_list as acc_f:
        for rec in SeqIO.parse(input_f, 'fasta'):

            if '>' in rec.description:
                headers = rec.description.split('>')
            else:
                headers = [rec.description]

            for header in headers:
                rec.id = extract_acc(header)
                rec.description = header
                out_f.write(rec.format('fasta'))
                acc_f.write(extract_acc(header) + "\n")


if __name__ == '__main__':
    main()
