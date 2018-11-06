#!/usr/bin/env python2

"""
Writen by Serafim Nenarokov.
6.11.2018

run: ./eukref_recover_headers.py --original <path to original fasta>
                                 --formatted <path to formatted fasta>
                                 -o <path to output>

Original file:
>ACC001 [Some genus species] - some metadata
GTAGTAGTAGTA
>ACC002 second metadata
CATCATCATCAT
>ACC003 third metadata
CATCATCATCAT

Formatted file:
>ACC001
GTAGTAGTAGTA
>ACC002
CATCATCATCAT
>OTHER_ACC042 something new here
GGGTACCCATATAAACA

The result of the script:
>ACC001 [Some genus species] - some metadata
GTAGTAGTAGTA
>ACC002 second metadata
CATCATCATCAT
>OTHER_ACC042 something new here
GGGTACCCATATAAACA
"""

import argparse

from Bio import SeqIO


def prepare_args():
    description = 'Recover the original fasta headers'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--original",
                        required=True,
                        type=argparse.FileType('r'),
                        help="Path to original fasta file.")

    parser.add_argument("--formatted",
                        required=True,
                        type=argparse.FileType('r'),
                        help="Path to formatted fasta file.")

    parser.add_argument("-o", "--output",
                        required=True,
                        type=argparse.FileType('w'),
                        help="Path where output fasta file will be saved.")

    return parser.parse_args()


def main():
    opts = prepare_args()

    orig_headers = {}

    with opts.original as orig_f:
        for rec in SeqIO.parse(orig_f, 'fasta'):
            orig_headers[rec.id] = rec.description

    with opts.formatted as format_f, opts.output as out_f:
        for rec in SeqIO.parse(format_f, 'fasta'):
            if rec.id in orig_headers:
                rec.description = orig_headers[rec.id]

            out_f.write(rec.format('fasta'))


if __name__ == '__main__':
    main()
