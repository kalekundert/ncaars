#!/usr/bin/env python3

import pyrosetta
import autoprop
import byoc
import logging

from pyrosetta import rosetta
from rosetta_utils import (
        DesignApp, Output1, init_rosetta, init_logging, xml_objects, kv
)
from byoc import Key, Method, DocoptConfig
from pathlib import Path

logger = logging.getLogger('design_fast_relax')

class DesignFastRelax(DesignApp, Output1):
    """\
Use the fast design algorithm to make an aaRS to bind a non-native NCAA.

Usage:
    design_fast_relax.py <pdb> <lig> [<scaffold>] [-o <pdb>] [-f <sfxn>]
        [-r <resfile>] [-p <pssm>] [-m <sele>] [-b <dist>] [-dD]

Arguments:
${app.shared_args}

Options:
${app.output_options}

${app.shared_options}

    -m --min-sele <sele>
        Which residues to minimize.  There are two options:

        shell:
            Only minimize residues that are either being repacked or are 
            adjacent (in primary sequence) to such a residue.  This is the 
            default.
        all:
            Minimize all residues.  This should be used in conjunction with the 
            `-b` option to prevent excess backbone movement.

    -b --bb-restraint <dist>
        Apply starting-coordinate restraints to the backbone atoms.  By 
        default, no restraints are applied.

Algorithm:
    FastRelax (also called FastDesign when mutations are allowed) is perhaps 
    the most common Rosetta design algorithm.  It works by iterating between 
    repacking and minimizing, all while slowly ramping score function weights 
    up and down.  Expect each design run to take about 1-2h.

References:
    [Nivón2013] DOI:10.1371/journal.pone.0059004
    [Conway2013] DOI:10.1002/pro.2389
    [Khatib2011] DOI:10.1073/pnas.1115898108
    [Tyka2011] DOI:10.1016/j.jmb.2010.11.008
"""

    min_sele = byoc.param(
            Key(DocoptConfig, '--min-sele'),
            default='shell',
    )
    bb_restraint = byoc.param(
            Key(DocoptConfig, '--bb-restraint'),
            default=None,
    )

    def main(self):
        super().main()

        rosetta_args = [
                '-in:file:extra_res_fa', self.lig_path,
        ]
        if self.bb_restraint is not None:
            rosetta_args += ['-relax:coord_cst_stdev', self.bb_restraint]

        init_logging(logger)
        init_rosetta(*rosetta_args, test_cycles=self.debug_run)

        pose = fast_relax(
                self.pose,
                scorefxn=self.scorefxn,
                resfile_path=self.resfile_path,
                pssm_path=self.pssm_path,
                min_sele=self.min_sele,
                bb_restraint=self.bb_restraint,
                dry_run=self.dry_run,
        )
        self.dump_pdb(pose, logger)

def fast_relax(pose, *, scorefxn, resfile_path, pssm_path, min_sele, bb_restraint, dry_run):
    xml_script = Path(__file__).parent / 'design_fast_relax.xml'
    script_vars = [
            f'scorefxn={scorefxn}',
            f'resfile_path={resfile_path}',
            f'pssm_path={pssm_path}',
            f'min=min_{min_sele}',
    ]
    xo = xml_objects(xml_script, script_vars)

    logger.info("Applying PSSM restraints to pose.")
    pssm = xo.get_mover('pssm')
    pssm.apply(pose)

    logger.info("Starting FastDesign protocol:")
    relax = xo.get_mover('relax')

    if dry_run:
        return pose

    if bb_restraint is not None:
        logger.info(kv("Restraints", "on"))
        logger.info(kv("Restrainted atoms", "backbone only"))
        logger.info(kv("Ramp restraints", "off"))
        logger.info(kv("Restraint stdev", bb_restraint))
        relax.constrain_relax_to_start_coords(True)
        relax.coord_constrain_sidechains(False)
        relax.ramp_down_constraints(False)
    else:
        logger.info(kv("Restraints", "off"))

    relax.apply(pose)
    return pose

if __name__ == '__main__':
    DesignFastRelax.entry_point()
