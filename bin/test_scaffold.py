#!/usr/bin/env python3

import parametrize_from_file
import prody

from rdkit import Chem
from scaffold import Scaffold, ConfigError
from voluptuous import Schema, Optional, Coerce, Or
from parametrize_from_file.voluptuous import Namespace
from pathlib import Path
from pprint import pprint
from pytest import approx, fixture

with_py = Namespace()
with_ncaars = Namespace(
        Scaffold=Scaffold,
        ConfigError=ConfigError,
        CWD=Path.cwd(),
        BUILTIN_SCAFFOLDS=Path(__file__).parent.parent / 'scaffolds',
)

@parametrize_from_file
def test_config_path(scaffold, expected):
    s = with_ncaars.exec(scaffold)['s']
    expected = with_ncaars.eval(expected)
    assert s.config_path == expected

@parametrize_from_file(
        schema=Schema({
            'scaffold': str,
            **with_ncaars.error_or({
                'expected': str,
            }),
        }),
)
def test_pdb_path(scaffold, expected, error):
    s = with_ncaars.exec(scaffold, get='s')
    expected = with_ncaars.eval(expected)

    with error:
        assert s.pdb_path == expected

@parametrize_from_file(
        schema=Schema({
            Optional('scaffold', default=''): str,
            Optional('files', default={}): {str: str},
            'expected': str,
        }),
        indirect=['files'],
)
def test_scorefxn(scaffold, files, expected):
    s = make_scaffold(scaffold, files)
    assert s.scorefxn == expected

@parametrize_from_file(
        schema=Schema({
            Optional('scaffold', default=''): str,
            Optional('files', default={}): {str: str},
            'expected_len': Coerce(int),
            'expected_apo_len': Coerce(int),
        }),
        indirect=['files'],
)
def test_atoms(scaffold, files, expected_len, expected_apo_len):
    s = make_scaffold(scaffold, files)
    assert len(s.atoms) == expected_len
    assert len(s.apo_atoms) == expected_apo_len

@parametrize_from_file(
        schema=Schema({
            Optional('scaffold', default=''): str,
            Optional('files', default={}): {str: str},
            Optional('pdb_ligands', default={}): Or({str: dict}, str),
            **with_ncaars.error_or({
                'expected': str,
            }),
        }),
        indirect=['files', 'pdb_ligands'],
)
def test_adenylate_mol_2d(scaffold, files, pdb_ligands, expected, error, monkeypatch):
    s = make_scaffold(scaffold, files)
    with error:
        assert Chem.MolToSmiles(s.adenylate_mol_2d) == expected

@parametrize_from_file(
        schema=Schema({
            Optional('scaffold', default=''): str,
            Optional('files', default={}): {str: str},
            Optional('pdb_ligands', default={}): Or({str: dict}, str),
            'expected': {
                'smiles': str,
                Optional('coords'): with_py.eval,
            },
        }),
        indirect=['files', 'pdb_ligands'],
)
def test_adenylate_mol_3d(scaffold, files, pdb_ligands, expected):
    s = make_scaffold(scaffold, files)

    assert Chem.MolToSmiles(s.adenylate_mol_3d) == expected['smiles']

    try:
        expected_coords = expected['coords']
    except KeyError:
        pass
    else:
        actual_coords = [
                tuple(v)
                for v in s.adenylate_mol_3d.GetConformer().GetPositions()
        ]

        pprint(expected_coords)
        pprint(actual_coords)
        assert s.adenylate_mol_3d.GetNumConformers() == 1
        assert actual_coords == approx(expected_coords)


@fixture
def files(request, tmp_path):
    for name, contents in request.param.items():
        p = tmp_path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(contents)

    return tmp_path

@fixture
def pdb_ligands(request, monkeypatch):
    def mock_fetch_pdb_ligand(resn):
        try:
            return request.param[resn]
        except KeyError:
            raise OSError

    if request.param != 'download':
        monkeypatch.setattr(prody, 'fetchPDBLigand', mock_fetch_pdb_ligand)

def make_scaffold(scaffold, files):
    if scaffold:
        return with_ncaars.exec(scaffold, get='s')
    else:
        return Scaffold(files)
