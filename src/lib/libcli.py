"""Argument parsing"""
import argparse

def parse_args():
    """
    Argument parser 
    """
    parser = argparse.ArgumentParser(
        description='Reduces the size of an MD \
        trajectory by removing water and selecting only CB and CA atoms.\
        All the other parameters are optional.')
    parser.add_argument('topology', type=str,
                        help='Topology file required for the trajectory')
    parser.add_argument('trajectory', type=str,
                        help='The trajectory')
    parser.add_argument('--mdtraj', action="store_true",
                        help='Force use of MDTraj')
    parser.add_argument('--mdanalysis', action="store_true",
                        help='Force use of MDAnalysis')
    parser.add_argument('--pytraj', action="store_true",
                        help='Force use of pytraj')
    parser.add_argument('--cpptraj', action="store_true",
                        help='Force use of CPPTRAJ')
    parser.add_argument('--gmx', action="store_true",
                        help='Force use of GROMACS')
    parser.add_argument('--vmd', action="store_true",
                        help='Force use of VMD')
    args = parser.parse_args()
    return args
