import os
import argparse
from eukref import utils


def format_arguments(args):
    args.output = os.path.normpath(args.output)
    args.db = os.path.normpath(args.db)
    return args


def parse_arguments():
    description = 'Genbank SSU rRNA gene parser'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input",
                        required=True,
                        type=argparse.FileType('r'),
                        help="Dataset of sarting sequences in fasta format.")

    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path do working directory")

    parser.add_argument("-d", "--db",
                        required=False,
                        default=utils.eukref_root('data', 'ncbi_nt'),
                        help="Path do folder with NCBI NT database")

    parser.add_argument("-r", "--ref",
                        required=True,
                        type=argparse.FileType('r'),
                        help="Path to Reference database, formatted as DNA, string")

    parsed = parser.parse_args()
    format_arguments(parsed)
    return parsed
