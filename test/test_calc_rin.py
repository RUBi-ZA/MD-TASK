#!/usr/bin/env python
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
    import calc_network
except ImportError or ModuleNotFoundError:
    print("ERROR: Module calc_network could not be found from src directory")
    sys.exit(1)

class Data:
    """MDM-TASK test data"""
    def __init__(self, topology_name, trajectory_name=None, stride=250):
        self.topology_name = topology_name
        self.trajectory_name = trajectory_name
        self.stride = stride
        if not trajectory_name:
            self.trajectory = md.load(self.topology_name)
        else:
            self.trajectory = md.load(self.trajectory_name, top=self.topology_name, stride=stride)
        self.args = calc_network.parse_args().parse_known_args()[0]

data = Data("../data/pr2.pdb", "../data/pr2.pdb")

class TestDRN(unittest.TestCase):
    """Test DRN metrics"""

    def get_centrality(self, metric):
        """
        Inputs:
         metric: [calc_BC, calc_CC, calc_DC, calc_EC, calc_ECC, calc_KC, calc_PR]
        """
        setattr(data.args, metric, True)
        value = calc_network.calc_centralities(
                    data.trajectory,
                    data.trajectory_name,
                    data.trajectory.n_frames, data.args)
        setattr(data.args, metric, False)
        return value

    def test_bc(self):
        """self explanatory"""
        val = self.get_centrality("calc_BC")['mean']["BC"].values[-1]
        self.assertGreater(val, -1)

    def test_cc(self):
        """self explanatory"""
        val = self.get_centrality("calc_CC")['mean']["CC"].values[-1]
        self.assertGreater(val, -1)

    def test_dc(self):
        """self explanatory"""
        val = self.get_centrality("calc_DC")['mean']["DC"].values[-1]
        self.assertGreater(val, -1)

    def test_ec(self):
        """self explanatory"""
        val = self.get_centrality("calc_EC")['mean']["EC"].values[-1]
        self.assertGreater(val, -1)

    def test_ecc(self):
        """self explanatory"""
        val = self.get_centrality("calc_ECC")['mean']["ECC"].values[-1]
        self.assertGreater(val, -1)

    def test_kc(self):
        """self explanatory"""
        val = self.get_centrality("calc_katz")['mean']["katz"].values[-1]
        self.assertGreater(val, -1)

    def test_pr(self):
        """self explanatory"""
        val = self.get_centrality("calc_PR")['mean']["PR"].values[-1]
        self.assertGreater(val, -1)

    def test_l(self):
        """self explanatory"""
        val = self.get_centrality("calc_L")['mean']["L"].values[-1]
        self.assertGreater(val, -1)

if __name__ == "__main__":
    unittest.main()
