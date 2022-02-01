#!/usr/bin/env python3

"""\
Copy an adenylated NCAA into an aaRS scaffold.

Usage:
    copy_ncaa_into_scaffold.py <scaffold> <ncaa> [-o <path>]

Options:
    -o --output <path>      [default: ncaa_adenylate_aars.pdb]
        The path where the model combined NCAA-aaRS model should be written.
"""

import prody
import docopt
from scaffold import Scaffold

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    scaffold = Scaffold(args['<scaffold>'])
    ncaa = prody.parsePDB(args['<ncaa>'])

    # Many protocols, including coupled moves, assume/require that the protein 
    # is chain A and the ligand is chain X.
    scaffold.apo_atoms.setChids('A')
    ncaa.setChids('X')

    # Give the ligand a defined residue number, so that it's easy to refer to 
    # in rosetta scripts.  In contrast, don't renumber the protein.  The 
    # scaffold includes config files (like resfiles) that may depend on the 
    # original numbering.
    ncaa.setResnums([1] * len(ncaa))

    out = scaffold.apo_atoms + ncaa
    prody.writePDB(args['--output'], out)

