#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.scaffold = ''
params.ncaa = ''
params.pdb = ''
params.sdf = ''
params.n_ncaa_confs = 10
params.out = 'results'

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
        path 'NCA.pdb', emit: pdb

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

/*
process design_greedyopt {
    input:
        path 'input.pdb'
        path 'NCA.params'

    output:
        path 'output.pdb'
}

process design_coupled_moves{
}
*/

workflow {
    scaffold = channel.value(params.scaffold)

    if( params.pdb ) {
        pdb = channel.fromPath(params.pdb)
    }
    else {
        if( params.sdf ) {
            sdf = channel.fromPath(params.sdf)
        }
        else {
            sdf = make_ncaa_conformers_etkdg(scaffold)
        }
        lig = molfile_to_params(sdf)
        pdb = copy_ncaa_into_scaffold(scaffold, lig.pdb)
    }
}
