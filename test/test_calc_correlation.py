#!/usr/bin/env python3
"""
calc_correlation testing
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
    import calc_correlation
except ImportError or ModuleNotFoundError:
    print("ERROR: Module ' calc_correlation.py' could not be found from src directory")
    sys.exit(1)

class Args():
    """MDM-TASK test data"""
    def __init__(self):
        self.topology = "../data/pr2.pdb"
        self.trajectory = "../data/pr2.xtc"
        self.step = 200
        self.select_atoms = "CA"
        self.lazy_load = False
        self.title = "Protein"
        self.prefix = "correlation"

class TestPRS(unittest.TestCase):
    """Test contact map"""

    def test_main(self):
        """Test whether contact map is of the right dimension"""
        var =  calc_correlation.main(Args())
        self.assertTupleEqual(var.shape, (198, 198))

if __name__ == "__main__":
    unittest.main()
