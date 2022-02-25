To make the starting model, I followed these steps:

- Load 1zh0 (the PDB-REDO version) and 5n5u (the original PDB version) into 
  pymol.

  - I decided to use the PDB-REDO version of 1zh0 because it's slightly better 
    in most metrics.  I don't think it really matters, though.  As far as I can 
    tell, both models are nearly identical other than for a few solvent exposed 
    rotamers.

  - I used the original PDB version of 5n5u because that's what I did for the 
    1j1u scaffold.  If I recall correctly, the PDB-REDO version of 5n5u wasn't 
    a clear improvement on the PDB version.

- Remove all water molecules:

    PyMol> remove resn HOH

- Remove the Tris.

    PyMol> remove resn TRS

- Remove all of 5n5u except the AMP.

    PyMol> remove 5n5u and not resn AMP

- Move the AMP and NpAla into chain X:

    PyMol> extract 1zh0_5n5u, all
    PyMol> alter resn AMP, chain='X'
    PyMol> alter resn AMP, resi=1
    PyMol> alter resn NAL, chain='X'
    PyMol> alter resn NAL, resi=2

- Save the resulting model to `1zh0_5n5u.pdb`.

- Relax the model using ``rosetta_relax_b``
