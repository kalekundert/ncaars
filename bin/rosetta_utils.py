#!/usr/bin/env python3

import pyrosetta
import logging

def init_rosetta(*flags, test_cycles=False, dry_run=False):
    flags = list(flags)

    if test_cycles:
        flags.append('-run:test_cycles')
    if dry_run:
        flags.append('-run:dry_run')

    pyrosetta.init(
            ' '.join(map(str, flags)),
            set_logging_handler='logging',
    )

def init_logging(logger):
    # The reason for using the logging module is to get uniform output 
    # formatting by taking advantage of the fact that pyrosetta can redirect 
    # all of rosetta's output to a logger.
    log = logging.getLogger()
    log.setLevel('WARNING')

    logger.setLevel('DEBUG')
    logging.getLogger('rosetta').setLevel('DEBUG')

    formatter = logging.Formatter('{asctime}\t{name}\t{levelname}\t{message}', style='{')
    handler = logging.StreamHandler()

    for handler in log.handlers[:]:
        log.removeHandler(handler)

    handler.setFormatter(formatter)
    log.addHandler(handler)

def kv(key, value):
    """Utility for logging key/value pairs in a consistent format."""
    return f"- {key+':':<35} {value}"


