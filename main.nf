#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.scaffold = ''  // If specified as a path, must be absolute.
params.ncaa = ''
params.pdb = ''
params.sdf = ''
params.out = 'results'
params.resfile = ''
params.n_ncaa_confs = 100
params.n_designs = 1000
params.dry_run = false

params.greedyopt = true
params.greedyopt_n = params.dry_run ? 2 : params.n_designs

params.fast_relax = true
params.fast_relax_n = params.dry_run ? 2 : params.n_designs

params.coupled_moves = true
params.coupled_moves_n = params.dry_run ? 2 : params.n_designs
params.coupled_moves_ligand_weight = 1
params.coupled_moves_bb_mover = 'backrub'

process make_ncaa_conformers_etkdg() {
    cpus 1
    memory '2GB'
    time '4h'
    label 'pyrosetta'
    publishDir params.out, mode: 'rellink'
    stageInMode 'rellink'

    input:
        val scaffold

    output:
        path 'ncaa_adenylate.sdf'

    """
    make_ncaa_conformers_etkdg.py \
        '${scaffold}' \
        '${params.ncaa}' \
        -n ${params.n_ncaa_confs} \
    """
}

process molfile_to_params {
    cpus 1
    memory '1GB'
    time '30m'
    label 'molfile_to_params'
    publishDir params.out, mode: 'rellink'
    stageInMode 'rellink'

    input:
        path 'ncaa_adenylate.sdf'

    output:
        tuple path('NCA.params'), path('NCA.pdb'), path('NCA_conformers.pdb')

    // References:
    // https://www.rosettacommons.org/demos/latest/tutorials/prepare_ligand/prepare_ligand_tutorial
    // https://www.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-ligands
    """
    molfile_to_params.py \
        -n NCA \
        -p NCA \
        --conformers-in-one-file \
        ncaa_adenylate.sdf
    """
}

process copy_ncaa_into_scaffold {
    cpus 1
    memory '1GB'
    time '30m'
    label 'pyrosetta'
    publishDir params.out, mode: 'rellink'
    stageInMode 'rellink'

    input:
        val scaffold
        tuple path('NCA.params'), path('NCA.pdb'), path('NCA_conformers.pdb')

    output:
        path 'ncaa_adenylate_aars.pdb'

    """
    copy_ncaa_into_scaffold.py '${scaffold}' NCA_conformers.pdb
    """
}

process design_fast_relax {
    cpus 1
    memory '2GB'
    time '8h'
    label 'pyrosetta'
    publishDir "${params.out}/fast_relax", mode: 'rellink'
    stageInMode 'rellink'

    input:
        path 'input.pdb'
        tuple path('NCA.params'), path('NCA.pdb'), path('NCA_conformers.pdb')
        val scaffold
        val n

    output:
        path "out_${n}.pdb"

    """
    design_fast_relax.py \
        input.pdb \
        NCA.params \
        '${scaffold}' \
        ${(r = params.resfile) ? "-r $r": ""} \
        ${params.dry_run ? "-D": ""} \
        -o 'out_${n}.pdb' \
    """
}

process design_coupled_moves {
    cpus 1
    memory '2GB'
    time '4h'
    label 'pyrosetta'
    publishDir "${params.out}/coupled_moves", mode: 'rellink'
    stageInMode 'rellink'

    input:
        path 'input.pdb'
        tuple path('NCA.params'), path('NCA.pdb'), path('NCA_conformers.pdb')
        val scaffold
        val n

    output:
        path "out_${n}.pdb"

    """
    design_coupled_moves.py \
        input.pdb \
        NCA.params \
        '${scaffold}' \
        ${(r = params.resfile) ? "-r $r": ""} \
        ${params.dry_run ? "-D": ""} \
        -l ${params.coupled_moves_ligand_weight} \
        -b ${params.coupled_moves_bb_mover} \
        -o 'out_${n}.pdb' \
    """
}

process design_greedyopt {
    cpus 1
    memory '2GB'
    time '48h'
    label 'pyrosetta'
    publishDir params.out, mode: 'rellink'
    stageInMode 'rellink'

    input:
        path 'input.pdb'
        tuple path('NCA.params'), path('NCA.pdb'), path('NCA_conformers.pdb')
        val scaffold

    output:
        path 'greedyopt_beyer2020'
        path 'GreedyOptTable.tab'

    """
    design_greedyopt_beyer2020.py \
        input.pdb \
        NCA.params \
        '${scaffold}' \
        ${(r = params.resfile) ? "-r $r": ""} \
        ${params.dry_run ? "-D": ""} \
        -n ${params.greedyopt_n} \
        -o 'greedyopt_beyer2020/{:03}.pdb.gz' \
    """
}

workflow {
    scaffold = channel.value(params.scaffold)

    if( params.pdb && file(params.pdb).exists() ) {
        pdb = channel.fromPath(params.pdb)
    }
    else {
        if( params.sdf && file(params.sdf).exists() ) {
            sdf = channel.fromPath(params.sdf)
        }
        else {
            sdf = make_ncaa_conformers_etkdg(scaffold)
        }
        lig = molfile_to_params(sdf) | first
        pdb = copy_ncaa_into_scaffold(scaffold, lig)
    }

    if( params.fast_relax ) {
        n = channel.of(1..params.fast_relax_n)
        design_fast_relax(pdb, lig, scaffold, n)
    }
    if( params.coupled_moves ) {
        n = channel.of(1..params.coupled_moves_n)
        design_coupled_moves(pdb, lig, scaffold, n)
    }
    if( params.greedyopt ) {
        design_greedyopt(pdb, lig, scaffold)
    }
}
