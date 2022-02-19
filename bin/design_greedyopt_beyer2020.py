#!/usr/bin/env python3

import pyrosetta
import autoprop
import byoc
import logging

from pyrosetta import rosetta, Vector1
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
from pyrosetta.rosetta.protocols.moves import MoverStatus
from rosetta_utils import DesignApp, OutputN, init_rosetta, init_logging, kv
from byoc import DocoptConfig
from pathlib import Path

logger = logging.getLogger('design_greedyopt_beyer2020')

@autoprop
class DesignGreedyOpt(DesignApp, OutputN):
    """\
Use the GreedyOpt algorithm described by [Beyer2020] to make an aaRS to bind
a non-native NCAA.

Usage:
    design_greedyopt_beyer2020.py <pdb> <lig> [<scaffold>] [-n <int>]
        [-o <pdb>] [-f <sfxn>] [-r <resfile>] [-dD]

Arguments:
${app.shared_args}

Options:
${app.output_options}

${app.shared_options}

Algorithm:
    This design protocol uses the Rosetta GreedyOpt algorithm [Nivon2014].  The 
    idea of using this algorithm and many of the specific parameters are taken 
    from [Beyer2020], which applies the GreedyOpt algorithm to the problem of 
    aaRS design.  However, this protocol is not meant to exactly recapitulate 
    the results from [Beyer2020]:

    - The original protocol was tailored to a specific aaRS, while this 
      one is general.
    - A different score function is used.
    - The repack shell is more conservative.

    It's also worth noting that GreedyOpt is meant to be used at the end of a 
    design simulation as an alternative to the "human intuition" that often 
    goes into choosing which designs to order, which mutations to revert, etc.  
    GreedyOpt is not meant to be a design simulation on its own.  However, 
    Beyer et al. were able to use it that way because they were starting with a 
    crystal structure of their aaRS already binding their NCAA of interest.  
    This structure came from a previous round of directed evolution, and in 
    effect the directed evolution served as the design simulation.

    The point of describing all this is to explain that there's reason to 
    expect that this algorithm will not perform well.  But it's still valuable 
    as a control.

References:
    [Beyer2020] DOI:10.1016/j.jmb.2020.06.014
    [Nivon2014] DOI:10.1002/prot.24463
"""
    default_n = 10

    def main(self):
        self.load(DocoptConfig)

        init_logging(logger)
        init_rosetta(
                '-extra_res_fa', self.lig_path,
                test_cycles=self.debug_run,
        )

        poses = iter_greedyopt_beyer2020(
                self.pose,
                scorefxn=self.scorefxn,
                resfile_path=self.resfile_path,
                dry_run=self.dry_run,
        )
        for i, pose in zip(range(self.n), poses):
            self.dump_pdb(pose, i, logger)

def iter_greedyopt_beyer2020(pose, *, scorefxn, resfile_path, dry_run):
    xml_script = Path(__file__).parent / 'design_greedyopt_beyer2020.xml'
    options = [
            f"scorefxn={scorefxn}",
            f"resfile_path={resfile_path}",
    ]

    rsp = RosettaScriptsParser()
    mover = rsp.generate_mover(str(xml_script), Vector1(options))

    if dry_run:
        Path('GreedyOptTable.tab').touch()
        while True:
            yield pose

    while True:
        pose_copy = Pose(pose)
        mover.apply(pose_copy)
        if mover.get_last_move_status() != MoverStatus.MS_SUCCESS:
            return
        yield pose_copy

if __name__ == '__main__':
    DesignGreedyOpt.entry_point()
