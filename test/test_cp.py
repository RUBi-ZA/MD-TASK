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
    import cp
except ImportError or ModuleNotFoundError:
    print("ERROR: Module 'cp.py' could not be found from src directory")
    sys.exit(1)

class Args():
    """MDM-TASK test data"""
    def __init__(self):
        self.topology = "../data/pr2.pdb"
        self.trajectory = "../data/pr2.xtc"
        self.step = 200
        self.num_frames = None
        self.prefix = "cp"

class TestCP(unittest.TestCase):
    """Test contact map"""

    def test_main(self):
        """Test whether contact map is of the right dimension"""
        val = cp.main(Args())
        self.assertTupleEqual(val.shape, (199, 199))

if __name__ == "__main__":
    unittest.main()
