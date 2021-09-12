#!/usr/bin/env python
""" Script for trajectory reduction to C-alphas and betas"""

__author__ = ["Olivier Sheik Amamuddy"]
__copyright__ = "Copyright 2020, Research Unit in Bioinformatics"
__license__ = "GNU GPL 3.0"
__version__ = "1.2"
__email__ = "oliserand@gmail.com"
__status__ = "Production"
__date__ = "22nd Nov 2020"

from lib.libcheck import check_tools
from lib.libcli import parse_args
from lib.libreduce import engines

def main(args):
    engine_selections = {"mdtraj":args.mdtraj, "MDAnalysis":args.mdanalysis,
                         "pytraj":args.pytraj, "cpptraj":args.cpptraj,
                         "gmx":args.gmx, "vmd":args.vmd}
    for engine_selected in engine_selections:
        if engine_selections[engine_selected]:
            print("INFO: Using {}.".format(engine_selected))
            engine = engine_selected
            break
    else:
        engine = check_tools()
    engines[engine](args.topology, args.trajectory)

if __name__ == "__main__":
    args = parse_args()
    main(args)
