#!/usr/bin/env python3

import pyrosetta
import autoprop
import byoc
import logging

from pyrosetta import rosetta, Vector1
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
from pyrosetta.rosetta.protocols.moves import MoverStatus
from rosetta_utils import init_rosetta, init_logging, kv
from scaffold import Scaffold
from byoc import Key, Method, DocoptConfig
from pathlib import Path

logger = logging.getLogger('design_greedyopt_beyer2020')

@autoprop
class DesignGreedyOpt(byoc.App):
    """\
Use the GreedyOpt algorithm described by [Beyer2020] to make an aaRS to bind
a non-native NCAA.

Usage:
    design_greedyopt_beyer2020.py <pdb> <lig> [<scaffold>] [-n <int>]
        [-o <pdb>] [-f <sfxn>] [-r <resfile>] [-d]

Arguments:
    <pdb>
        The starting model to design.  Chain A of this model must be an aaRS 
        scaffold, and chain X must be an NCAA-adenylate ligand.

    <lig>
        The parameter file for the NCAA-adenylate ligand.

    <scaffold>
        Either the name of a built-in scaffold, or the path to a custom 
        scaffold directory.  If specified, default values for various design 
        parameters will be taken from the scaffold. 

Options:
    -n <int>                    [default: 10]
        The number of structures to generate.

    -o --output <path>          [default: out_{:03}.pdb.gz]
        The path where the output files should be written.  The path should 
        contain python format code (e.g. '{}').  This will be used to add a 
        unique index to each path.

    -f --scorefxn <name>
        Which score function to use.  This overrides the score function 
        specified in the scaffold (i.e. the score function that the scaffold 
        was relaxed with), so this option should never be used for production 
        simulations.

    -r --resfile <path>
        The path to the resfile to use.  This overrides the default resfile 
        that may be specified in the scaffold.

    -d --dry-run
        Run an abbreviated simulation, for the purpose of debugging.  Note 
        however that this option doesn't do very much on its own (it just makes 
        the relax steps a bit faster).  To really abbreviate this protocol, it 
        is necessary to provide a resfile with very few (e.g. 2-5) allowed 
        mutations.

Algorithm:
    The design protocol uses the Rosetta GreedyOpt algorithm [Nivon2014].  The 
    idea of using this algorithm and many of the specific parameters are taken 
    from [Beyer2020], which applies the GreedyOpt algorithm to the problem of 
    aaRS design.  However, this algorithm is not meant to exactly recapitulate 
    the results from [Beyer2020]:

    - The original algorithm was tailored to a specific aaRS, while this 
      algorithm is general.
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
    __config__ = [DocoptConfig]

    pdb_path = byoc.param(
            Key(DocoptConfig, '<pdb>'),
    )
    lig_path = byoc.param(
            Key(DocoptConfig, '<lig>'),
    )
    n = byoc.param(
            Key(DocoptConfig, '-n'),
            cast=int,
    )
    scaffold = byoc.param(
            Key(DocoptConfig, '<scaffold>'),
            cast=Scaffold,
    )
    scorefxn = byoc.param(
            Key(DocoptConfig, '--scorefxn'),
            Method(lambda self: self.scaffold.scorefxn),
    )
    resfile_path = byoc.param(
            Key(DocoptConfig, '--resfile'),
            Method(lambda self: self.scaffold.resfile_path),
    )
    output_path = byoc.param(
            Key(DocoptConfig, '--output'),
    )
    dry_run = byoc.param(
            Key(DocoptConfig, '--dry-run'),
    )

    def main(self):
        self.load(DocoptConfig)

        init_rosetta('-extra_res_fa', self.lig_path, test_cycles=self.dry_run)
        init_logging(logger)

        results = iter_design_greedyopt_beyer2020(
                self.pose,
                scorefxn=self.scorefxn,
                resfile_path=self.resfile_path,
        )
        for i, pose in zip(range(self.n), results):
            out_path = self.output_path.format(i)
            Path(out_path).parent.mkdir(exist_ok=True, parents=True)
            pose.dump_pdb(out_path)

    def get_pose(self):
        return pose_from_pdb(self.pdb_path)


def iter_design_greedyopt_beyer2020(pose, *, scorefxn, resfile_path):
    xml_script = Path(__file__).parent / 'design_greedyopt_beyer2020.xml'
    options = [
            f"scorefxn={scorefxn}",
            f"resfile_path={resfile_path}",
    ]

    rsp = RosettaScriptsParser()
    mover = rsp.generate_mover(str(xml_script), Vector1(options))

    while True:
        pose_copy = Pose(pose)
        mover.apply(pose_copy)
        if mover.get_last_move_status() != MoverStatus.MS_SUCCESS:
            return
        yield pose_copy

if __name__ == '__main__':
    DesignGreedyOpt.entry_point()
