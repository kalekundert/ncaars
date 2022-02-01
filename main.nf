#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.scaffold = ''
params.ncaa = ''
params.pdb = ''
params.sdf = ''
params.out = 'results'
params.dry_run = false
params.n_ncaa_confs = 10
params.resfile = ''

params.beyer = true
params.beyer_n = 10

process make_ncaa_conformers_etkdg() {
    conda "${baseDir}/conda/rdkit_prody.yml"
    publishDir params.out, mode: 'rellink'

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
    container 'ncaars-molfile-to-params:v2021.49-dev61812'
    publishDir params.out, mode: 'rellink'

    input:
        path 'ncaa_adenylate.sdf'

    output:
        path 'NCA.params', emit: param
        path 'NCA_conformers.pdb', emit: pdb

    // References:
    // https://www.rosettacommons.org/demos/latest/tutorials/prepare_ligand/prepare_ligand_tutorial
    // https://www.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-ligands
    """
    molfile_to_params.py \
        -n NCA \
        -p NCA \
        --conformers-in-one-file \
        ncaa_adenylate.sdf

    # Don't know why `molfile_to_params.py` doesn't do this for me...
    mv NCA.pdb NCA_conformers.pdb
    """
}

process copy_ncaa_into_scaffold {
    conda "${baseDir}/conda/rdkit_prody.yml"
    publishDir params.out, mode: 'rellink'

    input:
        val scaffold
        path 'NCA.pdb'

    output:
        path 'ncaa_adenylate_aars.pdb'

    """
    copy_ncaa_into_scaffold.py '${scaffold}' NCA.pdb
    """
}

process design_greedyopt {
    container 'ncaars-rdkit-prody-pyrosetta'
    publishDir params.out, mode: 'rellink'

    input:
        path 'input.pdb'
        path 'NCA.params'
        path 'NCA_conformers.pdb'
        val scaffold

    output:
        path 'design_greedyopt_beyer2020'
        path 'GreedyOptTable.tab'

    """
    design_greedyopt_beyer2020.py \
        input.pdb \
        NCA.params \
        '${scaffold}' \
        ${(r = params.resfile) ? "-r $r": ""} \
        ${params.dry_run ? "-d": ""} \
        -n ${params.beyer_n} \
        -o 'design_greedyopt_beyer2020/{:03}.pdb.gz' \
    """
}

/*
process design_coupled_moves{
}
*/

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
        lig = molfile_to_params(sdf)
        pdb = copy_ncaa_into_scaffold(scaffold, lig.pdb) | first
    }

    if( params.beyer ) {
        design_greedyopt(pdb, lig, scaffold)
    }
}
