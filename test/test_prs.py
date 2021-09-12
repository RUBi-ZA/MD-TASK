#!/usr/bin/env python3
"""
DRN testing
"""
import sys
import unittest

try:
    import mdtraj as md
except ImportError or ModuleNotFoundError:
    print("ERROR: MDTraj module not found")
    sys.exit(1)
try:
    sys.path.append("../src")
    import prs
except ImportError or ModuleNotFoundError:
    print("ERROR: Module 'prs.py' could not be found from src directory")
    sys.exit(1)

class Args():
    """MDM-TASK test data"""
    def __init__(self):
        self.topology = "../data/pr2.pdb"
        self.trajectory = "../data/pr2.xtc"
        self.final = "../data/pr7.pdb"
        self.step = 200
        self.perturbations = 10
        self.num_frames = None
        self.aln = False
        self.prefix = "result"

class TestPRS(unittest.TestCase):
    """Test contact map"""

    def test_main(self):
        """Test whether contact map is of the right dimension"""
        shape = prs.main(Args())
        self.assertTupleEqual(shape, (198, ))

if __name__ == "__main__":
    unittest.main()
