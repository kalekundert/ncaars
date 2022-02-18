To make the starting model, I followed these steps:

- Load 1j1u (the PDB-REDO version) and 5n5u (the original PDB version) into 
  pymol.

- Superimpose 5n5u on 1j1u.

- Remove all water molecules.

- Remove the magnesium.

- Remove the tRNA (1j1u_redo and chain B).

- Remove all of 5n5u except the AMP.

- Save the resulting model to `1j1u_5n5u.pdb`.

- Manually edit the resulting file:

  - Put Y401 and AMP in chain X.
  - Terminate chain A after L306.

- Renumber the atoms consecutively.

- Relax the model using ``rosetta_relax_b``
