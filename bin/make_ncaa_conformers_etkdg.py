#!/usr/bin/env python3

"""\
Generate conformers for an adenylated non-canonical amino acid (NCAA) using the 
ETKDG algorithm (via RDKit).

Usage:
    make_adenylated_ncaa_conformers.py <scaffold> <ncaa> [-o <sdf>]
        [-n <confs>] [-N <tries>] [-t <rms>] [-a <tol>] [-p <tol>] [-k <frac>]
        [-r <seed>] [-i <ca>] 

Arguments:
    <scaffold>
        Either the name of a built-in scaffold (e.g. 'mma_pylrs', 'mja_tyrrs') 
        or a path to a directory containing all the necessary information about 
        a custom scaffold.

    <ncaa>
        A SMILES string describing the NCAA.

Options:
    -o --output <sdf>                   [default: ncaa_adenylate.sdf]
        The path where the SDF file containing the generated conformers will be 
        written.

    -n --num-confs <int>                [default: 10]
        The number of conformations to generate.

    -N --num-tries-per-conf <int>       [default: -1]
        Most attempts to generate a conformation are unsuccessful, most often  
        because not all the atoms are close enough to their required positions.  
        This option controls the number of attempts to make before giving up.  
        Specifically, the total number of iterations attempted will be the 
        product of `--num-confs` and `--num-tries-per-conf`.  Once this limit 
        is exceeded, the program will terminate with however many conformations 
        have been found to that point.  By default, the program will not 
        terminate until it has found the requested number of conformations.

    -t --rms-threshold <float>          [default: -1]
        Every conformer must be at least this distance from every other 
        conformer.  The default value (-1) indicates that no threshold will be 
        applied.

    -a --anchor-tol <float>             [default: 0.5]
        How close each anchor atom must be to its target position.  Any 
        conformations with greater distances than this will be discarded.

    -p --pocket-tol <float>             [default: 3.0]
        How close each NCAA atom must be to a corresponding pocket atom defined 
        in the scaffold.  Pocket atoms are meant to loosely identify where the 
        binding pocket is located.  Correspondence between NCAA and pocket 
        atoms are defined based on the number of bonds that separate those 
        atoms from a common reference atom (e.g. Cα).  It's possible for each 
        NCAA atom to correspond with multiple pocket atoms, in which case it 
        only needs to be within the given distance for one of them.

    -k --quantile-keep <frac>           [default: 0.9]
        Keep only the highest-scoring fractions of all the generated 
        conformations that pass all the distance cutoffs.  The purpose of this 
        is to avoid the possibility that a handful of outlier conformations 
        with really bad geometries get incorporated into the rotamer set.

    -r --random-seed <seed>             [default: 1]
        The seed used to randomize atom coordinates during the distance 
        geometry step of the algorithm.

    -i --nucleophile <index>
        The index of the NCAA atom to attach AMP to.  If the NCAA contains a 
        single α-amino acid moiety, this atom will be identified automatically.  
        If the NCAA either contains multiple α-amino acid moieties or has a 
        non-standard backbone, this index needs to be specified.  Note that 
        atoms are numbered in the order they appear in the SMILES string (i.e. 
        going from left to right), starting at 0 and skipping any hydrogens.

References:
    [Riniker2015] DOI:10.1021/acs.jcim.5b00654
"""

import sys
import tidyexc
import numpy as np

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.DistanceGeometry.DistGeom import DoTriangleSmoothing
from scaffold import mol_from_smiles, mol_from_smarts, smiles_from_mol, smarts_from_mol
from dataclasses import dataclass
from itertools import product, combinations, count
from more_itertools import first, zip_equal
from log import init_logging, log_call, log_message
from copy import copy
from math import ceil

class SubstrateAlignment:

    @classmethod
    def from_scaffold(cls, mol, scaffold):
        return cls.from_substructs(
                mol,
                scaffold.adenylate_mol_3d, 
                scaffold.adenylate_anchor_query,
                scaffold.adenylate_pocket_query,
        )

    @classmethod
    def from_substructs(cls, mol, ref, anchor_substruct, pocket_substruct):
        mol_names = {
                mol: 'adenylated NCAA',
                ref: 'scaffold',
        }
        substruct_names = {
                anchor_substruct: 'anchor',
                pocket_substruct: 'pocket',
        }
        def require_substruct(mol, substruct):
            indices = list(mol.GetSubstructMatch(substruct))

            if not indices:
                err = UsageError(
                        mol=mol,
                        mol_name=mol_names[mol],
                        substruct=substruct,
                        substruct_name=substruct_names[substruct],
                )
                err.brief = "can't find {substruct_name} atoms in {mol_name}"
                err.info += lambda e: f"{e.mol_name}: {smiles_from_mol(e.mol)}"
                err.info += lambda e: f"{e.substruct_name}: {smarts_from_mol(e.substruct)}"
                raise err

            return indices

        mol_anchor_indices = require_substruct(mol, anchor_substruct)
        ref_anchor_indices = require_substruct(ref, anchor_substruct)

        anchor_indices = dict(zip(
                mol_anchor_indices,
                ref_anchor_indices,
        ))
        ref_pocket_indices = require_substruct(ref, pocket_substruct)

        return cls(mol, ref, anchor_indices, ref_pocket_indices)

    def __init__(self, mol, ref, anchor_indices, ref_pocket_indices):
        self.mol = mol
        self.ref = ref

        conformers = ref.GetConformers()
        if len(conformers) != 1:
            err = UsageError(conformers=conformers)
            err.brief = lambda e: f"scaffold adenylate must have only 1 conformer, found: {len(e.conformers)}"
            raise err
        self.ref_conf = first(conformers)

        self.anchor_indices = anchor_indices
        self.anchor_coords = {
                i: self.ref_conf.GetAtomPosition(j)
                for i, j in anchor_indices.items()
        }
        self.pocket_indices = [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetSymbol() != 'H' and \
                   atom.GetIdx() not in anchor_indices
        ]

        self.ref_anchor_indices = {
                v: k for k, v in self.anchor_indices.items()
        }
        self.ref_pocket_indices = {}

        topo_dists_mol = Chem.GetDistanceMatrix(self.mol)
        topo_dists_ref = Chem.GetDistanceMatrix(self.ref)

        # Pick the anchor atom nearest to the pocket atoms.  This is necessary 
        # to handle the case where some of the anchor atoms are on a different 
        # molecule than the pocket atoms.

        min_topo_dist = np.inf

        for k_ref in self.ref_anchor_indices:
            dists = topo_dists_ref[k_ref, list(ref_pocket_indices)]
            dist = np.sum(dists)

            if dist < min_topo_dist:
                i_ref = k_ref
                min_topo_dist = dist

        i_mol = self.ref_anchor_indices[i_ref]

        for j_ref in ref_pocket_indices:
            self.ref_pocket_indices[j_ref] = set()
            d = topo_dists_ref[i_ref, j_ref]

            for j_mol in self.pocket_indices:
                if topo_dists_mol[i_mol, j_mol] == d:
                    self.ref_pocket_indices[j_ref].add(j_mol)

        assert any(self.ref_pocket_indices.values())



class UsageError(tidyexc.Error):
    pass

class SamplingError(tidyexc.Error):
    pass

@log_call
def find_nucleophile(mol):
    bb = mol_from_smarts('OC(=O)CN')
    hits = mol.GetSubstructMatches(bb)

    if len(hits) == 0:
        err = UsageError(mol=mol, bb=bb)
        err.brief = "can't find α-amino acid backbone"
        err.info += lambda e: f"NCAA: {smiles_from_mol(e.mol)}"
        err.info += lambda e: f"query: {smarts_from_mol(e.bb)}"
        err.hints += "use the '--nucleophile <index>' option to specify where AMP should be attached"
        raise err

    if len(hits) > 1:
        err = UsageError(mol=mol, bb=bb)
        err.brief = "found multiple α-amino acid moieties"
        err.info += lambda e: f"NCAA: {smiles_from_mol(e.mol)}"
        err.info += lambda e: f"query: {smarts_from_mol(e.bb)}"
        err.info += lambda e: f"AMP could be attached to the following atom indices: {', '.join(repr(h[0]) for h in hits)}"
        err.hints += "use the '--nucleophile <index>' option to specify where AMP should be attached"
        raise err

    return hits[0][0]

@log_call
def attach_amp(ncaa, nuc_i):
    # From wikipedia: https://en.wikipedia.org/wiki/Adenosine_monophosphate
    amp_smiles = 'P(O)(=O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3O'
    amp = mol_from_smiles(amp_smiles)
    n = len(ncaa.GetAtoms())
    m = len(amp.GetAtoms())

    adenylate = Chem.RWMol()
    adenylate.InsertMol(ncaa)
    adenylate.InsertMol(amp)
    adenylate.AddBond(nuc_i, n, Chem.BondType.SINGLE)

    Chem.SanitizeMol(adenylate)

    for i in range(n):
        adenylate.SetProp('ncaarsType', 'ncaa')
    for i in range(n, n + m):
        adenylate.SetProp('ncaarsType', 'amp')

    return adenylate.GetMol()

@log_call
def generate_conformers_etkdg(
        alignment, *,
        num_confs=1,
        random_seed=0,
        rms_thresh=-1,
        anchor_restraint_tol=0.01,
        anchor_prune_tol=0.5,
        pocket_prune_tol=3.0,
        max_tries_per_conf=0,
        quantile_keep=0.9,
):
    mol_h = Chem.AddHs(alignment.mol)

    # Try making my own distance matrix
    bounds = AllChem.GetMoleculeBoundsMatrix(mol_h)
    DoTriangleSmoothing(bounds)
    restrain_anchor_atoms(alignment, bounds, tol=anchor_restraint_tol)

    params = AllChem.srETKDGv3()
    params.SetBoundsMat(bounds)
    params.clearConfs = False
    params.randomSeed = random_seed

    # Don't use the `pruneRMSThresh` option.  It includes the restrained atoms 
    # in its calculation, and so gets rid of way too many conformations.

    n_build = int(ceil(num_confs / quantile_keep))
    n_failed_embedding = 0
    n_rejected_anchor = 0
    n_rejected_pocket = 0
    n_rejected_score = 0
    n_kept = 0

    scores = {}

    if max_tries_per_conf > 0:
        counter = range(n_build * max_tries_per_conf)
    else:
        counter = count(0)

    for i in counter:
        params.randomSeed += 1
        conf_id = AllChem.EmbedMolecule(mol_h, params)
        if conf_id < 0:
            n_failed_embedding += 1
            continue

        conf = mol_h.GetConformer(conf_id)

        AllChem.AlignMol(
                mol_h, alignment.ref, conf_id,
                atomMap=list(alignment.anchor_indices.items()),
        )

        if any_anchor_atoms_misplaced(alignment, conf, tol=anchor_prune_tol):
            mol_h.RemoveConformer(conf_id)
            n_rejected_anchor += 1
            continue

        if any_pocket_atoms_misplaced(alignment, conf, tol=pocket_prune_tol):
            mol_h.RemoveConformer(conf_id)
            n_rejected_pocket += 1
            continue

        sfxn = AllChem.MMFFGetMoleculeForceField(
                mol_h,
                AllChem.MMFFGetMoleculeProperties(mol_h),
                confId=conf_id,
        )
        scores[conf_id] = sfxn.CalcEnergy()

        n_kept += 1
        if n_kept >= n_build:
            break

    if scores:
        score_cutoff = np.quantile(
                list(scores.values()),
                quantile_keep,
        )
        for conf_id, score in scores.items():
            if score > score_cutoff:
                mol_h.RemoveConformer(conf_id)
                n_rejected_score += 1
                n_kept -= 1

    log_message(f"attempted to generate {i+1} conformers")
    log_message(f"failed to generate a conformer {n_failed_embedding} times")
    log_message(f"discarded {n_rejected_anchor} conformers due to misplaced anchor atoms")
    log_message(f"discarded {n_rejected_pocket} conformers due to misplaced pocket atoms")
    log_message(f"discarded the {n_rejected_score} worst-scoring conformers")
    log_message(f"kept {n_kept} conformers")

    if n_kept == 0:
        err = SamplingError(mol=alignment.mol)
        err.brief = "failed to generate conformers"
        err.info += lambda e: f"mol: {smiles_from_mol(e.mol)}"
        err.hints += "see log for details"
        raise err

    return mol_h

@log_call
def restrain_anchor_atoms(alignment, bounds, *, tol=0.01):
    anchor_coords = alignment.anchor_coords

    for i, j in combinations(sorted(anchor_coords), 2):
        # `i` is guaranteed to be lower than `j`
        vi = anchor_coords[i]
        vj = anchor_coords[j]
        d = vi.Distance(vj)
        bounds[i,j] = d + tol
        bounds[j,i] = max(d - tol, 0)

def any_anchor_atoms_misplaced(alignment, conf, *, tol=0.5):
    if tol < 0:
        return False

    for i, xyz_expected in alignment.anchor_coords.items():
        xyz_actual = conf.GetAtomPosition(i)
        if xyz_expected.Distance(xyz_actual) > tol:
            return True

    return False

def any_pocket_atoms_misplaced(alignment, conf, *, tol=3.0):
    if tol < 0:
        return False

    for i in alignment.ref_pocket_indices:
        xyz_expected = alignment.ref_conf.GetAtomPosition(i)

        for j in alignment.ref_pocket_indices[i]:
            xyz_actual = conf.GetAtomPosition(j)
            if xyz_expected.Distance(xyz_actual) < tol:
                break
        else:
            return True

    return False

@log_call
def write_sdf(mol, path):
    with Chem.SDWriter(path) as sdf:
        for conformer in mol.GetConformers():
            sdf.write(mol, conformer.GetId())
    

if __name__ == '__main__':
    import docopt
    from scaffold import Scaffold

    init_logging()
    args = docopt.docopt(__doc__)

    scaffold = Scaffold(args['<scaffold>'])
    ncaa = mol_from_smiles(smiles := args['<ncaa>'])

    if not ncaa:
        err = UsageError(smiles=smiles)
        err.brief = "can't parse NCAA SMILES string: {smiles}"
        raise err

    if x := args['--nucleophile']:
        nuc_i = int(x)
    else:
        nuc_i = find_nucleophile(ncaa)

    ncaa_adenylate_2d = attach_amp(ncaa, nuc_i)
    ncaa_adenylate_alignment = SubstrateAlignment.from_scaffold(
            ncaa_adenylate_2d,
            scaffold,
    )
    ncaa_adenylate_3d = generate_conformers_etkdg(
            ncaa_adenylate_alignment,
            num_confs=int(args['--num-confs']),
            random_seed=int(args['--random-seed']),
            rms_thresh=float(args['--rms-threshold']),
            anchor_prune_tol=float(args['--anchor-tol']),
            pocket_prune_tol=float(args['--pocket-tol']),
            max_tries_per_conf=int(args['--num-tries-per-conf']),
            quantile_keep=float(args['--quantile-keep']),
    )
    write_sdf(ncaa_adenylate_3d, args['--output'])

