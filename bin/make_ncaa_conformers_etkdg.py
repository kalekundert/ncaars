#!/usr/bin/env python3

"""\
Generate conformers for an adenylated non-canonical amino acid (NCAA) using the 
ETKDG algorithm (via RDKit).

Usage:
    make_adenylated_ncaa_conformers.py <scaffold> <ncaa> [-o <sdf>] [-i <ca>] 
        [-n <confs>] [-t <rms>] [-r <seed>]

Arguments:
    <scaffold>
        Either the name a built-in scaffolds (e.g. 'mma_pylrs', 'mja_tyrrs') or 
        a path to a directory containing all the necessary information about a 
        custom scaffold.

    <ncaa>
        A SMILES string describing the NCAA.

Options:
    -o --output <sdf>               [default: ncaa_adenylate.sdf]
        The path where the SDF file containing the generated conformers will be 
        written.

    -n --num-confs <int>            [default: 10]
        The number of conformations to generate.

    -t --rms-threshold <float>      [default: -1]
        Every conformer must be at least this distance from every other 
        conformer.  The default value (-1) indicates that no threshold will be 
        applied.

    -r --random-seed <seed>         [default: 1]
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

# TODO:
# - Pocket restraints.  I disabled them because they seems to distort bond 
#   geometries the way I have them written now.  I think I need to do something 
#   where I generate models freely, then throw them out if they don't look like 
#   how I want them to.  That will let me do comparisons in absolute space, 
#   instead of pairwise-distance-space.

import sys
import tidyexc
import pandas as pd

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.DistanceGeometry.DistGeom import DoTriangleSmoothing
from scaffold import mol_from_smiles, mol_from_smarts, smiles_from_mol, smarts_from_mol
from dataclasses import dataclass
from itertools import product, combinations
from more_itertools import zip_equal
from log import init_logging, log_call, log_message
from copy import copy

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
        self.ref_conf = next(iter(conformers))

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
        self.ref_pocket_indices = ref_pocket_indices

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
        pocket_restraint_tol=1.0,
        anchor_prune_tol=0.5,
):
    mol_h = Chem.AddHs(alignment.mol)

    # Try making my own distance matrix
    bounds = AllChem.GetMoleculeBoundsMatrix(mol_h)
    DoTriangleSmoothing(bounds)
    restrain_anchor_atoms(alignment, bounds)
    #restrain_pocket_atoms(alignment, bounds)

    params = AllChem.srETKDGv3()
    params.SetBoundsMat(bounds)
    params.randomSeed = random_seed

    if random_seed == 0:
        err = SamplingError("random seed set to 0")
        err.info = "with the random seed is set to 0, every generated conformation would be identical"
        raise err

    # Don't use the `pruneRMSThresh` option.  It includes the restrained atoms 
    # in its calculation, and so gets rid of way too many conformations.

    conf_ids = AllChem.EmbedMultipleConfs(mol_h, num_confs, params)

    if not conf_ids:
        err = SamplingError(mol=mol)
        err.brief = "failed to generate conformers"
        err.info += lambda e: f"mol: {smailes_from_mol(e.mol)}"
        raise err

    log_message(f"generated {len(conf_ids)} conformers")

    # From: rdkit/Code/GraphMol/DistGeomHelpers/Embedder.h
    # coordMap    a map of int to Point3D, between atom IDs and their locations
    #             their locations.  If this container is provided, the
    #             coordinates are used to set distance constraints on the
    #             embedding. The resulting conformer(s) should have distances
    #             between the specified atoms that reproduce those between the
    #             points in \c coordMap. Because the embedding produces a
    #             molecule in an arbitrary reference frame, an alignment step
    #             is required to actually reproduce the provided coordinates.

    AllChem.AlignMol(
            mol_h, alignment.ref,
            atomMap=list(alignment.anchor_indices.items()),
    )
    AllChem.AlignMolConformers(
            mol_h,
            atomIds=list(alignment.anchor_indices.keys()),
    )

    n_pruned = prune_misaligned_conformers(mol_h, alignment)
    n_kept = len(mol_h.GetConformers())

    log_message(f"discarded {n_pruned} mis-aligned conformers")
    log_message(f"kept {n_kept} conformers")

    if n_kept == 0:
        err = SamplingError(mol=mol, n_pruned=n_pruned)
        err.brief = "all {n_pruned} generated conformers were pruned due to poor restraint satisfaction"
        err.info += lambda e: f"mol: {smailes_from_mol(e.mol)}"
        err.hint += "this may indicate the wrong chirality was specified for one or more atoms"
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

@log_call
def restrain_pocket_atoms(alignment, bounds, *, tol=1.0):
    pocket_dists = calc_pocket_dists(alignment)
    topo_dists = Chem.GetDistanceMatrix(alignment.mol)

    def upper(i, j):
        # Upper bounds go in the upper-right half of the bounds matrix.
        return tuple(sorted([i, j]))

    def lower(i, j):
        # Lower bounds go in the lower-left half of the bounds matrix.
        return upper(i, j)[::-1]

    for i in alignment.pocket_indices:
        for j in alignment.anchor_indices:
            d = topo_dists[i,j]
            p = pocket_dists.query('mol_anchor_i == @j and dist_topo == @d')

            if p.empty:
                continue

            bounds[upper(i,j)] = max(p.dist_3d) + tol
            bounds[lower(i,j)] = min(p.dist_3d) - tol


    return bounds

@log_call
def calc_pocket_dists(alignment):
    topo_dists_ref = Chem.GetDistanceMatrix(alignment.ref)
    conf = alignment.ref_conf
    dists = []

    for (i_mol, i_ref), j_ref in product(
            alignment.anchor_indices.items(),
            alignment.ref_pocket_indices,
    ):
        pi = conf.GetAtomPosition(i_ref)
        pj = conf.GetAtomPosition(j_ref)
        dists.append({
            'mol_anchor_i': i_mol,
            'ref_anchor_i': i_ref,
            'ref_pocket_i': j_ref,
            'dist_topo': topo_dists_ref[i_ref, j_ref],
            'dist_3d': pi.Distance(pj),
        })

    return pd.DataFrame(dists)

@log_call
def prune_misaligned_conformers(mol, alignment, *, tol=0.5):
    n_pruned = 0
    for conf in list(mol.GetConformers()):
        for i, xyz_expected in alignment.anchor_coords.items():
            xyz_actual = conf.GetAtomPosition(i)
            if xyz_expected.Distance(xyz_actual) > tol:
                n_pruned += 1
                mol.RemoveConformer(conf.GetId())
                break
    return n_pruned

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
    )
    write_sdf(ncaa_adenylate_3d, args['--output'])

