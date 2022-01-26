#!/usr/bin/env python3

import make_ncaa_conformers_etkdg as mnce
import pytest, parametrize_from_file
import numpy as np

from rdkit import Chem, Geometry
from scaffold import mol_from_smiles, mol_from_smarts
from dataclasses import asdict
from voluptuous import Schema, Optional, Coerce
from parametrize_from_file import star
from parametrize_from_file.voluptuous import Namespace
from pytest import approx
from pytest_unordered import unordered

with_py = Namespace()
with_mnce = Namespace(star(mnce))
with_math = Namespace('from math import *')

class MockSubstrateAlignment:
    pass

def alignment_from_mols(structs):
    return mnce.SubstrateAlignment.from_substructs(
            mol_from_smiles(structs['mol']),
            Chem.MolFromMolBlock(structs['ref']),
            mol_from_smarts(structs['anchor']),
            mol_from_smarts(structs['pocket']),
    )

def mock_alignment_from_attrs(attrs):
    alignment = MockSubstrateAlignment()

    try:
        alignment.mol = mol_from_smiles(attrs['mol'])
    except KeyError:
        pass

    try:
        alignment.anchor_indices = attrs['anchor_indices']
    except KeyError:
        pass

    try:
        alignment.anchor_coords = points_from_tuples(attrs['anchor_coords'])
    except KeyError:
        pass

    try:
        alignment.pocket_indices = attrs['pocket_indices']
    except KeyError:
        pass

    try:
        alignment.ref_pocket_indices = attrs['ref_pocket_indices']
    except KeyError:
        pass

    return alignment

alignment_from_mols.schema = {
        'mol': str,
        'ref': str,
        'anchor': str,
        'pocket': str,
}
mock_alignment_from_attrs.index_schema = {
        'anchor_indices': {Coerce(int): Coerce(int)},
        'anchor_coords': {Coerce(int): with_py.eval},
        'pocket_indices': [Coerce(int)],
        'ref_pocket_indices': [Coerce(int)],
}
mock_alignment_from_attrs.schema = {
        Optional('mol'): str,
        **{
            Optional(k): v
            for k, v in mock_alignment_from_attrs.index_schema.items()
        },
}

@parametrize_from_file(
        schema=Schema({
            'ncaa': str,
            **with_mnce.error_or({
                'expected': Coerce(int),
            }),
        }),
)
def test_find_nucleophile(ncaa, expected, error):
    ncaa = mol_from_smiles(ncaa)
    with error:
        assert mnce.find_nucleophile(ncaa) == expected

def test_attach_amp():
    l_ala = Chem.MolFromSmiles('C[C@@H](C(=O)O)N')
    adenylate = mnce.attach_amp(l_ala, 4)

    # The Cα stereochemistry indicator is switched because I swapped the order 
    # of the amine/carboxylate groups.  I confirmed manually (i.e. by looking 
    # at 2D structures in jupyter) that this expected SMILES string is correct.
    expected = Chem.MolFromSmiles('C[C@H](N)C(=O)OP(O)(=O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3O')

    assert Chem.MolToSmiles(adenylate) == Chem.MolToSmiles(expected)

@pytest.mark.parametrize('num_confs', [2, 20])
@parametrize_from_file(
        schema=Schema({
            'alignment': alignment_from_mols.schema,
            'positions': {Coerce(int): with_py.eval},
            'chiral_centers': {Coerce(int): str},
        }),
)
def test_generate_conformers_etkdg(alignment, num_confs, positions, chiral_centers):
    alignment = alignment_from_mols(alignment)

    def tuple_from_point(v):
        return v.x, v.y, v.z

    mol = mnce.generate_conformers_etkdg(
        alignment,
        num_confs=num_confs,
        random_seed=1,
    )
    assert len(mol.GetConformers()) == num_confs

    # Check positioning:

    for i, expected in positions.items():
        expected = Geometry.Point3D(*expected)

        for conf in mol.GetConformers():
            actual = conf.GetAtomPosition(i)

            # The positioning isn't perfect, so make sure that each fixed 
            # conformer atom is within 0.2Å of where it should be.
            assert actual.Distance(expected) < 0.2

    # Check stereochemistry:

    for conf in mol.GetConformers():
        Chem.AssignStereochemistryFrom3D(mol, conf.GetId())
        assert dict(Chem.FindMolChiralCenters(mol)) == chiral_centers

    # Check random seed:

    mol_same = mnce.generate_conformers_etkdg(
        alignment,
        num_confs=num_confs,
        random_seed=1,
    )
    mol_diff = mnce.generate_conformers_etkdg(
        alignment,
        num_confs=num_confs,
        random_seed=2,
    )
    for i in range(num_confs):
        for j in alignment.pocket_indices:
            p = mol.GetConformer(i).GetAtomPosition(j)
            ps = mol_same.GetConformer(i).GetAtomPosition(j)
            pd = mol_diff.GetConformer(i).GetAtomPosition(j)
            assert tuple_from_point(p) == approx(tuple_from_point(ps))
            assert tuple_from_point(p) != approx(tuple_from_point(pd))

@parametrize_from_file(
        schema=Schema({
            'alignment': alignment_from_mols.schema,
            **with_mnce.error_or({
                'expected': mock_alignment_from_attrs.index_schema,
            })
        }),
)
def test_substrate_alignment(alignment, expected, error):
    with error:
        a = alignment_from_mols(alignment)

    if not error:
        t = tuples_from_points
        x = mock_alignment_from_attrs(expected)

        assert a.anchor_indices == x.anchor_indices
        assert t(a.anchor_coords) == approx(t(x.anchor_coords))
        assert a.pocket_indices == x.pocket_indices
        assert a.ref_pocket_indices == x.ref_pocket_indices

@parametrize_from_file(
        schema=Schema({
            'alignment': mock_alignment_from_attrs.schema,
            'bounds': np.genfromtxt,
            'tol': Coerce(float),
            'expected': np.genfromtxt,
        }),
)
def test_restrain_anchor_atoms(alignment, bounds, tol, expected):
    alignment = mock_alignment_from_attrs(alignment)
    mnce.restrain_anchor_atoms(alignment, bounds, tol=tol)
    assert bounds == approx(expected)

@parametrize_from_file(
        schema=Schema({
            'alignment': alignment_from_mols.schema,
            'bounds': np.genfromtxt,
            'tol': Coerce(float),
            'expected': np.genfromtxt,
        }),
)
def test_restrain_pocket_atoms(alignment, bounds, tol, expected):
    alignment = alignment_from_mols(alignment)
    mnce.restrain_pocket_atoms(alignment, bounds, tol=tol)
    assert bounds == approx(expected)

@parametrize_from_file(
        schema=Schema({
            'alignment': alignment_from_mols.schema,
            'expected': with_math.eval,
        }),
)
def test_calc_pocket_dists(alignment, expected):
    alignment = alignment_from_mols(alignment)
    df = mnce.calc_pocket_dists(alignment)
    expected = unordered(approx(x) for x in expected)
    assert df.to_dict('records') == expected

def points_from_tuples(d):
    return {k: Geometry.Point3D(*v) for k, v in d.items()}

def tuples_from_points(d):
    return {k: (p.x, p.y, p.z) for k, p in d.items()}
