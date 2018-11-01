#!/usr/bin/env python3.6
"""eukref_gbretrieve.py main script.

Authors: Serafim Nenarokov, Martin Kolisko.
"""

import logging

from eukref import settings
from eukref.session import Session

TRIM_MIN_LEN = 500
MAX_ITERS = 256


def prepare_logger():
    """Initialize logger."""
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=logging.INFO)
    global logger
    logger = logging.getLogger()


def main():
    """Execute gbretrieve logic."""
    prepare_logger()

    args = settings.parse_arguments()
    settings.format_arguments()

    session = Session(session_path=args.output)
    session.add_starting_set(args.input)

    starting_set = session.starting_set()

    # running formatter & ucheme & trim short
    starting_set.normalize_headers()
    starting_set.uchime_and_save(args.ref)
    starting_set.filter_short(TRIM_MIN_LEN)
    starting_set.save()

    # main loop
    for loop_id in range(MAX_ITERS):
        logger.info(f"Performing iteration {loop_id}")

        iteration = session.init_next_iteration()
        iteration.add_dataset(starting_set.path)

        # do smth with this dataset


if __name__ == '__main__':
    main()
