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
    import contact_map
except ImportError or ModuleNotFoundError:
    print("ERROR: Module 'contact_map.py' could not be found from src directory")
    sys.exit(1)

class Data:
    """MDM-TASK test data"""
    def __init__(self, topology, trajectory_name=None, stride=250, residue="ASH25"):
        self.topology = topology
        self.trajectory_name = trajectory_name
        self.stride = stride
        self.residue = residue
        self.args = contact_map.parse_args()
        self.args.residue = self.residue
        self.args.topology = self.topology
        self.args.trajectory = self.trajectory_name

data = Data(topology="../data/pr2.pdb", trajectory_name="../data/pr2.xtc")

class TestContactMap(unittest.TestCase):
    """Test contact map"""

    def test_get_contact_map(self):
        """Test whether contact map is of the right dimension"""
        df = contact_map.get_contact_map(data.args)
        self.assertTupleEqual(df.shape, (12,3))

    def test_get_chain_labels(self):
        """Test natoms and unique chains"""
        val = contact_map.get_chain_labels(data.topology)
        natoms = len(val)
        nchains = len(set(val))
        self.assertTupleEqual((nchains, natoms), (3, 3194))

if __name__ == "__main__":
    unittest.main()
