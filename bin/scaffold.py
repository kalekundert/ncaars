#!/usr/bin/env python3

import autoprop
import byoc
import tidyexc
import prody

from rdkit.Chem import AllChem as Chem
from byoc import NtConfig
from more_itertools import unique_everseen as unique
from pathlib import Path
from io import StringIO

NAMED_SCAFFOLD_DIR = Path(__file__).parent.parent / 'scaffolds'

# Naming conventions:
# *_atoms: A `prody.AtomGroup` object
# *_mol: An `rdkit.Mol` object
# *_query: An `rdkit.Mol` object derived from a SMARTS string

@autoprop.cache
class Scaffold:
    """
    Provide access to the parameters defining a scaffold for the aaRS design 
    pipeline.

    The same scaffold can be designed to bind many different non-canonical 
    amino acids (NCAAs), so these parameters are deliberately agnostic to the 
    NCAA itself.  Instead, they specify:

    - A 3D aaRS structure containing a bound adenylate.
    - A selection (using the ProDy syntax) identifying the adenylate.
    - The 2D chemical structure of the adenylate.
    - Which adenylate atoms to hold fixed during the design simulations.
    - Which score function was used to relax the aaRS structure.

    Note that the adenylate in the structure would ideally have the native 
    amino acid attached, although it can be any amino acid (or none).

    This information is contained in a directory with the following structure:

    - ``scaffold.pdb``: 3D coordinates
    - ``config.nt``: all other information
    """
    __config__ = [
            NtConfig.setup(
                path_getter=lambda self: self.config_path,
            ),
    ]

    pdb_relpath = byoc.param(
            ('scaffold', 'pdb'),
            default='scaffold.pdb',
    )
    scorefxn = byoc.param(
            ('scaffold', 'scorefxn'),
    )
    resfile_relpath = byoc.param(
            ('scaffold', 'resfile'),
            default='resfile',
    )
    pssm_relpath = byoc.param(
            ('scaffold', 'pssm'),
            default='pssm',
    )
    frag_relpaths = byoc.param(
            ('scaffold', 'frags'),
            cast=lambda d: {int(k): v for k, v in d.items()}
    )
    adenylate_sele = byoc.param(
            ('adenylate', 'sele'),
    )
    adenylate_smiles = byoc.param(
            ('adenylate', 'smiles'),
            default=None,
    )
    adenylate_inchi = byoc.param(
            ('adenylate', 'inchi'),
            default=None,
    )
    adenylate_anchor_smarts = byoc.param(
            ('adenylate', 'anchor_smarts'),
    )
    adenylate_pocket_smarts = byoc.param(
            ('adenylate', 'pocket_smarts'),
    )

    def __init__(self, path_or_name):
        self._root = NAMED_SCAFFOLD_DIR / path_or_name
        if not self._root.is_dir() or not self._root.parent == NAMED_SCAFFOLD_DIR:
            self._root = Path(path_or_name).resolve()

    def __repr__(self):
        return f'{self.__class__.__name__}({str(self._root)!r})'

    def get_config_path(self):
        return self._root / 'config.nt'

    def get_pdb_path(self):
        pdb_path = self._root / self.pdb_relpath

        if pdb_path.suffix != '.pdb':
            err = ConfigError(self, pdb_path=pdb_path)
            err.brief = "PDB path must have '.pdb' extension"
            err.info += "PDB path: {pdb_path}"
            raise err

        return pdb_path

    def get_resfile_path(self):
        return self._root / self.resfile_relpath

    def get_pssm_path(self):
        return self._root / self.pssm_relpath

    def get_frag_paths(self):
        return {
                k: self._root / v
                for k, v in self.frag_relpaths.items()
        }

    def get_frag_flags(self):
        sizes = []
        paths = []

        # Sort the fragments by decreasing size of the fragments, because
        # rosetta insists that the fragment arguments be in this order.

        by_size = lambda x: -x[0]
        for size, path in sorted(self.frag_paths.items(), key=by_size):
            sizes.append(size)
            paths.append(path)

        # If no size-1 fragments were generated, but larger fragments were,
        # also add the 'none' pseudo-path.  This will cause rosetta to make
        # size-1 fragments from the next largest fragment set.

        if sizes and sizes[-1] > 1:
            sizes.append(1)
            paths.append('none')

        # Construct the command-line flags.

        flags = []

        if paths and sizes:
            flags.append('-loops:frag_sizes')
            flags.extend(map(str, sizes))
            flags.append('-loops:frag_files')
            flags.extend(map(str, paths))

        return flags

    def get_atoms(self):
        return prody.parsePDB(str(self.pdb_path))

    def get_apo_atoms(self):
        sele = self.adenylate_sele
        return self.atoms.select(f'not ({sele})').toAtomGroup()

    def get_adenylate_atoms(self):
        sele = self.adenylate_sele
        atoms_sele = self.atoms.select(sele)

        if atoms_sele is None:
            err = ConfigError(self, sele=sele)
            err.brief = "adenylate selection is empty"
            err.info += "PDB path: {scaffold.pdb_path}"
            err.info += "selection: {sele}"
            raise err

        return atoms_sele.toAtomGroup()

    def get_adenylate_mol_2d(self):
        if self.adenylate_smiles and self.adenylate_inchi:
            raise ConfigError(self, "adenylate specified twice (SMILES and InChI)")

        if x := self.adenylate_smiles:
            mol = mol_from_smiles(x.strip())
            if not mol:
                raise ConfigError(self, "can't parse SMILES string: {smiles}", smiles=x)
            return mol

        if x := self.adenylate_inchi:
            mol = Chem.MolFromInchi(x.strip())
            if not mol:
                raise ConfigError(self, "can't parse InChI string: {inchi}", inchi=x)
            return mol

        res = self.adenylate_atoms
        resns = list(unique(res.getResnames()))

        with ConfigError.add_info(
                "PDB path: {scaffold.pdb_path}"
                "selection: {sele}",
                sele=self.adenylate_sele,
                scaffold=self,
        ):
            if len(resns) > 1:
                err = ConfigError(self, resns=resns)
                err.brief = "adenylate selection matches multiple residues"
                err.info += lambda e: f"matching residues: {', '.join(e.resns)}"
                raise err

            resn = resns.pop()

            try:
                info = prody.fetchPDBLigand(resn)
            except OSError:
                err = ConfigError(self, resn=resn)
                err.brief = "unknown adenylate residue"
                err.info += "queried PDB ligand database for {resn!r}: not found"
                raise err

            try:
                inchi = info['InChI_InChI']
            except KeyError:
                err = ConfigError(self, resn=resn)
                err.brief = "unknown adenylate 2D chemical structure"
                err.info += "queried PDB ligand database for {resn!r}:\nresidue: found\nInChI string: not found"
                raise err

            mol = Chem.MolFromInchi(inchi.strip())
            if not mol:
                err = ConfigError(self, resn=resn, inchi=inchi)
                err.brief = "can't parse InChI string: {inchi}"
                err.info += "downloaded the above string from PDB ligand database for {resn!r}"
                raise err

            return mol

    def get_adenylate_mol_3d(self):
        pdb = StringIO()
        prody.writePDBStream(pdb, self.adenylate_atoms)
        mol = Chem.MolFromPDBBlock(pdb.getvalue())

        Chem.AssignStereochemistryFrom3D(mol)

        # This function generates a warning because there are two oxygens 
        # bonded to the phosphate and it's ambiguous which is which.  
        mol = Chem.AssignBondOrdersFromTemplate(self.adenylate_mol_2d, mol)

        return mol

    def get_adenylate_anchor_query(self):
        return mol_from_smarts(self.adenylate_anchor_smarts)

    def get_adenylate_pocket_query(self):
        return mol_from_smarts(self.adenylate_pocket_smarts)


class ConfigError(tidyexc.Error):

    def __init__(self, scaffold, *args, **kwargs):
        super().__init__(*args, **{'scaffold': scaffold, **kwargs})
        self.info += "config: {scaffold.config_path}"

def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol.SetProp('user_smiles', smiles, True)
    return mol

def smiles_from_mol(mol):
    try:
        return mol.GetProp('user_smiles')
    except KeyError:
        return Chem.MolToSmiles(mol)

def mol_from_smarts(smarts):
    mol = Chem.MolFromSmarts(smarts)
    if mol:
        mol.SetProp('user_smarts', smarts, True)
    return mol

def smarts_from_mol(mol):
    try:
        return mol.GetProp('user_smarts')
    except KeyError:
        return Chem.MolToSmarts(mol)


