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

def renumber_residues(atoms):
    resi = atoms.getResindices() + 1
    atoms.setResnums(resi)

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    scaffold = Scaffold(args['<scaffold>'])
    ncaa = prody.parsePDB(args['<ncaa>'])

    # Many protocols, including coupled moves, assume/require that the protein 
    # is chain A and the ligand is chain X.
    scaffold.apo_atoms.setChids('A')
    ncaa.setChids('X')

    out = scaffold.apo_atoms + ncaa
    renumber_residues(out)
    prody.writePDB(args['--output'], out)

