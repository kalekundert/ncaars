******
NCaaRS
******

.. image:: https://img.shields.io/badge/nextflow-DSL2-informational
   :alt: Nextflow version

.. image:: https://img.shields.io/readthedocs/ncaars.svg
   :alt: Documentation
   :target: https://ncaars.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/github/workflow/status/kalekundert/ncaars/Test%20and%20release/master
   :alt: Test status
   :target: https://github.com/kalekundert/ncaars/actions

.. image:: https://img.shields.io/github/last-commit/kalekundert/ncaars?logo=github
   :alt: Last commit
   :target: https://github.com/kalekundert/ncaars

Notes
=====

Installing the pipeline
-----------------------
- Either docker or conda are required.

  - Docker:

    - You must build the docker images yourself, because the build process 
      requires access to the rosetta source code and the pyrosetta conda 
      channel (both proprietary).

    - Scripts to build the images are found in the docker directory.

  - Conda:

    - I don't know why, but Nextflow seems unable to install the pyrosetta 
      conda package.  It does nothing for the amount of time specified by 
      `conda.createTimeout`, prints a stack traces, then continues to hang 
      indefinitely.
    - If you want to use conda (e.g. if docker is unavailable), you can get 
      around this by creating your own environments using the files in the 
      conda directory, and configuring nextflow to use those.

Running the pipeline
--------------------
- To see logs from conformer-generation steps:

    > conda install eliot-tree
    > eliot-tree log.json

Custom scaffolds
----------------
- Have to do this by hand.
- The basic steps:

  - Get a PDB model of the scaffold with a ligand bound.
  - Relax the model in a rosetta score function, e.g. ref2015.

    - I use kalekundert/rosetta_relax_b for this purpose.
    - I recommend renumbering the residues in the scaffold to count from 1 
      before this step.

  - Describe the adenylate ligand in the config file.

Custom design algorithms
------------------------
- Implement a `--dry-run` option; it's very useful for development.

