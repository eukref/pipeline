#!/usr/bin/env python3.6

import os
import sys

from eukref import settings
from eukref.session import Session
from eukref.eukref_fasta_file import EukrefFastaFile


def main():
    args = settings.parse_arguments()

    session = Session(session_path=args.output)

    # formatting
    starting_set = EukrefFastaFile(args.input)
    starting_set.normalize_headers()
    starting_set.save(session.starting_set_path())






if __name__ == '__main__':
    main()
