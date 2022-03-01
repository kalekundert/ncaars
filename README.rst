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
- Install nextflow

  - This is preinstalled on many clusters.

- Either docker or conda are required.  Docker is preferred:

  - Docker:

    You must build the docker images yourself, because the build process 
    requires access to the rosetta source code and the pyrosetta conda channel 
    (both proprietary).  Here's an outline of the process:

    - Get access to the Rosetta source code:

      - If you are a developer, just configure git to be able clone the main 
        Rosetta repository without needing to provide any passwords (e.g.  
        upload your SSH key to GitHub)

      - If you are not a developer: get a Rosetta license, download the source 
        code, and move it into the ``docker/molfile_to_params/rosetta`` 
        directory.  The following path should exist::
          
          docker/molfile_to_params/rosetta/source/scripts/python/public/molfile_to_params.py

    - Get access to the PyRosetta conda channel::
        
        https://<user>:<pass>@conda.graylab.jhu.edu

    - Build all images from scratch::

        $ cd docker
        $ PYROSETTA_CONDA_CHANNEL=... make

    - Many cluster systems require singularity containers, instead of docker 
      containers.  You can create these as follows::

        $ make sif

      - Note that the ``ncaars_rdkit_prody_pyrosetta`` image is very big (6 
        GB).  If your ``/tmp`` is not big enough, you may need to specify a 
        different temporary directory like so::

          SINGULARITY_TMPDIR=/home/kale make sif

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
- Make a ``nextflow.config`` file:

  - Specify which NCAA you want to design for::

      params.ncaa = '...'  // SMILES

  - Specify which aaRS scaffold you want to use::

      params.scaffold = 'mma_pylrs'

    There are two built-in scaffolds you can use: ``mma_pylrs`` and 
    ``mja_tyrrs``.  The following section describes how to create a custom 
    scaffold.
      
  - Tell nextflow which containers/conda environments to use:

    - Each process is labeled either "pyrosetta" or "molfile_to_params", 
      indicating which environment is needed.

    - Use these labels to provide the appropriate container names for your      
      installation::

        process {
          withLabel: pyrosetta {
            container = 'ncaars-rdkit-prody-pyrosetta:latest'
          }
          withLabel: molfile_to_params {
            container = 'ncaars-molfile-to-params:v2021.49-dev61812'
          }
        }

- Run nextflow::

    $ nextflow kalekundert/main.nf

  There is no need to actually download this repository to run the pipeline.  
  The above command will instruct nextflow to automatically download and run 
  the pipeline.  However, if you have the repository downloaded anyways (e.g.  
  to build the docker images), you can also give the path to the ``main.nf`` 
  script.

- The ``make_ncaa_conformers_etkdg.py`` script generates the following 
  warning::

    WARNING: More than one matching pattern found - picking one

  This happens because the AMP phosphate has two identical oxygens, so RDKit 
  can't unambiguously match the oxygens in the SMARTS strings to the oxygens in 
  the PDB file.  As long as no effort is made to treat the two oxygens 
  differently, though, it doesn't matter how they are matched.

- To see logs from conformer-generation steps:

    > conda install eliot-tree
    > eliot-tree log.json

Custom scaffolds
----------------
Preparing a custom scaffold for this pipeline takes a lot of setup work, and 
most of it has to be done by hand.  Below is an outline of the basic steps:

- Get a PDB model of the aaRS scaffold with the amino acid adenylate bound.  It 
  may be necessary to merge coordinates from two models to make this happen.  
  The exact coordinates of the AMP moiety are the most important; the amino 
  acid coordinates are just used to loosely define the binding pocket.

- Specify which ligand atoms to consider part of the "anchor" and the "pocket".  
  Anchor atoms are held in place throughout the design process, and are meant 
  to correspond to the AMP.  Pocket atoms are meant to loosely define where the 
  active site is, and are typically the sidechain atoms of the natural amino 
  acid substrate.

  Both sets of atoms are specified using SMARTS queries.  These queries can be 
  finicky and hard to get right.  I write them using a Jupyter lab session, 
  because that makes it easy to experiment (in no small part because rdkit 
  automatically shows 2D molecular structures in Jupyter sessions).  Here's an 
  example session::

    > cd /path/to/ncaars/bin
    > from scaffold import Scaffold
    > s = Scaffold('my_custom_scaffold')
    > m = s.adenylate_mol_2d
    > m
    2D structure of adenylate
    > from rdkit.Chem import AllChem as Chem
    > smarts = Chem.MolFromSmarts
    > q = smarts('OC(=O)')
    > q
    2D structure of query
    > m.GetSubstructMatch(q)
    list of matching positions
    > m
    2D structure of adenylate, this time with matching atoms highlighted

- Specify which residues will be allowed to mutate.  This is done using a 
  "resfile"; a rosetta-specific file format described here: 

  https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles
  
  It's recommended that you do not limit which amino acids are allowed at the 
  positions you want to design, and that you don't freeze any positions.  The 
  design algorithms will make these decisions themselves (e.g. limiting amino 
  acids based on secondary structure or a PSSM, freezing residues based on 
  their proximity to the design shell) and it's best not to step on their toes.

  For example, here is a resfile that allows design at positions 32 and 34::

    NATAA
    START

    32 A ALLAA
    34 A ALLAA
  
- Relax the model in a rosetta score function, e.g. ref2015.

  - I use kalekundert/rosetta_relax_b for this purpose.
  - I recommend renumbering the residues in the scaffold to count from 1 
    before this step.

- Create a FASTA file:

  - I did this by loading the PDB into PyMOL and using the ``save`` command 
    to make a FASTA file.

  - This file isn't directly used by the design pipeline, but it's needed to 
    make some of the other input files.

- Create a PSSM:

  - Used by design algorithms to bias towards stable sequences.

  - Don't provide an automatic script for this, because it requires the BLAST 
    database (specifically nr).  This is far too big to include in a docker 
    container, and unnecessary since most institutions already make the BLAST 
    databases available somehow.  So I'll just give the command here::

      psiblast \
          -db nr_v5 \
          -query 2zim.fasta \
          -out_ascii_pssm 2zim.pssm \
          -num_iterations 4 \
          -num_alignments 1 \
          -num_threads 8 \

    I request the following resources when running this command:

    - CPUs: 8
    - Memory: 80 GB (10 GB/core)
    - Time: 12h

- Create a fragment library:

  - Create an account on: https://old.robetta.org
  - Submit a "Fragment Library" job.
  - Upload the FASTA file created above.
  - Don't exclude homologues.  That option is only used for benchmarking.
  - You can use ``contrib/wget_robetta.sh`` to download the results.
      
Custom design algorithms
------------------------
- Most design algorithms take at least these arguments:

  - The path to a PDB model of the scaffold with the target NCAA in the binding 
    site.  The model will have been relaxed in the Rosetta force field in the 
    context of its native ligand.  The native ligand will have been replaced 
    with the target ligand without any further optimization, so there may be 
    severe clashes.  It is assumed that these clashes will be resolved by the 
    design algorithm itself.

  - The path to the Rosetta ligand parameter file for the target NCAA.  This 
    file should be provided to Rosetta via the ``-extra_res_fa`` command line 
    option.

  - The path to (or name of) the scaffold.  The scaffold contains a number of 
    default design parameters described in the section above.

  - `--dry-run` and `--debug-run` options; they're very useful for development.

