#!/usr/bin/env python3

from eliot.json import EliotJSONEncoder
from eliot import to_file, log_call, log_message

class ReprJsonEncoder(EliotJSONEncoder):
    SUPPORTED_TYPES = str, int, float, dict, list, tuple

    def default(self, obj):
        # Import within this function to avoid circular imports
        from rdkit import Chem

        if isinstance(obj, Chem.Mol):
            return Chem.MolToSmiles(obj)

        try:
            obj = super().default(obj)
        except TypeError:
            return repr(obj)


def init_logging(path='log.json', mode='w'):
    to_file(open('log.json', 'w'), encoder=ReprJsonEncoder)

