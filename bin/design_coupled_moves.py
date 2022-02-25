#!/usr/bin/env python3

import pyrosetta
import pandas as pd
import autoprop
import byoc
import logging
import re

from pyrosetta import rosetta, Vector1
from rosetta_utils import (
        DesignApp, Output1, init_rosetta, init_logging, xml_objects, kv,
)
from byoc import Key, Method, DocoptConfig
from pathlib import Path
from dataclasses import dataclass

logger = logging.getLogger('design_coupled_moves')

@autoprop
class DesignCoupledMoves(DesignApp, Output1):
    """\
Use the coupled moves algorithm to make an aaRS to bind a non-native NCAA.

Usage:
    design_coupled_moves.py <pdb> <lig> [<scaffold>] [-n <int>] [-o <pdb>]
        [-f <sfxn>] [-r <resfile>] [-p <pssm>] [-b <mover>] [-l <weight>] [-dD]

Arguments:
${app.shared_args}

Options:
${app.output_options}

${app.shared_options}

    -b --backbone-mover <name>
        Which algorithm to use when making backbone moves.  The following 
        options are recognized:

        backrub
            Small rotations based on movements seen in high-resolution crystal 
            structure.

        kic/fragment
            Copy backbone torsions from homologous fragments, and use KIC to 
            keep the bond lengths and angles correct.  This option requires 
            that the scaffold includes a fragment library.

        kic/walking
            Perturb backbone torsions by small amounts, and use KIC to keep the 
            all bond lengths and angles correct.

    -l --ligand-weight
        The score function weight for protein-ligand interactions.  Typically 
        set between 1.0-2.0, see [Ollikainen2015].

Algorithm:
    This design protocol uses the Rosetta coupled moves algorithm 
    [Ollikainen2015].  Compared to the other design algorithms implemented as 
    part of this pipeline, coupled moves puts particular emphasis on modeling 
    backbone movement.

    The central idea of this algorithm is to propose changes to both the 
    backbone and the sidechains before deciding whether to accept or reject the 
    proposed changes.  This should allow bigger sidechain moves to be made.  
    Redesigning the specificity of ligand binding sites was the problem that 
    coupled moves was originally designed for and benchmarked on 
    [Ollikainen2015], so there's reason to believe that it will perform well 
    for aaRS design.

References:
    [Ollikainen2015] DOI:10.1371/journal.pcbi.1004335
"""

    backbone_mover = byoc.param(
            Key(DocoptConfig, '--backbone-mover'),
            default='backrub',
    )
    ligand_weight = byoc.param(
            Key(DocoptConfig, '--ligand-weight', cast=float),
            default=1.0,
    )

    def main(self):
        self.load(DocoptConfig)

        # A lot of coupled-moves options can only be set via the command line.
        init_logging(logger)
        init_rosetta(
                '-extra_res_fa', self.lig_path,
                '-coupled_moves::ligand_mode', True,
                '-coupled_moves::number_ligands', 1,
                '-coupled_moves::ligand_weight', self.ligand_weight,
                '-out::mute', 'protocols.backrub.BackrubMover',
                *self.backbone_mover_flags,
        )

        # Add a logging handler that will let us record each sequence visited 
        # by coupled moves.  This is a bit hacky, but unfortunately there's no 
        # other way to get this information out of rosetta.

        traj_handler = TrajectoryHandler()
        logging.getLogger('rosetta').addHandler(traj_handler)

        # This is described more in my lab notebook, but:
        #
        # - It probably takes more than 10,000 iterations for mutations at each 
        #   position to be sampled.
        #
        # - After ≈20,000 iterations, scores start *increasing* dramatically.
        #
        # Based on this, I decided to set a hard limit of 10,000 iterations.  I 
        # think that's about the most I can do before running into problems.

        pose = coupled_moves(
                self.pose,
                scorefxn=self.scorefxn,
                resfile_path=self.resfile_path,
                pssm_path=self.pssm_path,
                n_trials=10000 if not self.debug_run else 100,
                dry_run=self.dry_run,
        )
        self.dump_pdb(pose, logger)

    def get_backbone_mover_flags(self):
        if self.backbone_mover == 'backrub':
            return [
                '-coupled_moves::backbone_mover', 'backrub',
            ]

        if self.backbone_mover == 'kic/fragment':
            return [
                '-coupled_moves::backbone_mover', 'kic',
                '-coupled_moves::kic_perturber', 'fragment',
                *self.scaffold.frag_flags,
            ]

        if self.backbone_mover == 'kic/walking':
            return [
                '-coupled_moves::backbone_mover', 'kic',
                '-coupled_moves::kic_perturber', 'walking',
            ]

        raise ValueError(f"unknown backbone mover: {self.backbone_mover}")

class TrajectoryHandler(logging.Handler):

    @dataclass
    class Frame:
        i: int
        seq: str
        score: float

    def __init__(self):
        super().__init__()
        self.trajectory = []

    @property
    def trajectory_df(self):
        return pd.DataFrame(self.trajectory)

    def emit(self, record):
        # Note: any exceptions that occur in this method will silently ignored 
        # by the logging framework!

        msg = record.getMessage()
        pattern = r'protocols.CoupledMovesProtocol: (\x1b\[0m)?(\d+) ([A-Z]+) (-?\d+\.\d+)'

        m = re.search(pattern, msg)
        if m := re.search(pattern, msg):
            frame = self.Frame(
                    i=int(m.group(2)),
                    seq=m.group(3),
                    score=float(m.group(4)),
            )
            self.trajectory.append(frame)

        # I can also wait until I see:
        #   Design Positions: 108 111 112 115 150 152 154 190 205 223 225 227 Starting Sequence: ALMAQCAYSAGG
        # 
        # Then I know how long the sequence should be, and can make the pattern 
        # more specific.  Not sure it's worth it, though.
        #
        # pattern = r'protocols.CoupledMovesProtocol: (\d)+ ([A-Z]+) (-?\d+\.\d+)'
        # pattern = 'protocols.CoupledMovesProtocol: 998 ALMAQQAHQAGG -470.465'

    def createLock(self):
        self.lock = None

def coupled_moves(pose, *, scorefxn, resfile_path, pssm_path, n_trials, dry_run):
    xml_script = Path(__file__).parent / 'design_coupled_moves.xml'
    script_vars = [
            f"scorefxn={scorefxn}",
            f"resfile_path={resfile_path}",
            f"pssm_path={pssm_path}",
    ]
    xo = xml_objects(xml_script, script_vars)

    if dry_run:
        return pose

    pssm = xo.get_mover('pssm')
    pssm.apply(pose)

    cm = xo.get_mover('coupled_moves')
    cm.set_ntrials(n_trials)
    cm.apply(pose)

    # This is the lowest scoring pose, not the last pose.
    return pose

if __name__ == '__main__':
    DesignCoupledMoves.entry_point()
