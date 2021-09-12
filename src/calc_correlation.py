#!/usr/bin/env python
"""
Calculates Dynamic Cross Correlations from MD data
"""
#
# Calculate correlation in MD trajectory
# Script distributed under GNU GPL 3.0
# Author: Caroline Ross

__version__ = 1.2
__date__ = "15th November 2020"

import sys
import math
import argparse
import numpy as np
import matplotlib
from lib.cli import CLI
from lib.utils import Logger
from lib.trajectory import load_trajectory
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_traj(traj, topology=None, step=1, selected_atoms="CA",
               lazy_load=False):
    """Reads a trajectory and returns a dictionary of residues with their
    coordinates from each frame reshaped as nested lists"""
    selected_atoms =  selected_atoms.strip().split(",")
    traj = load_trajectory(traj, topology, step, lazy_load)[0]
    residues = {}
    for frame in traj:
        for atom in frame.topology.atoms:
            if atom.name in selected_atoms:
                res = atom.residue.index
                ac = frame.xyz[0, atom.index]
                co_ords = [ac[0], ac[1], ac[2]]
                if res in residues:
                    residues[res].append(co_ords)
                else:
                    residues[res] = [co_ords]
    return residues

def mean_dot(m1, m2, size):
    """Computes the pairwise dot product between two arrays of deviations
    for all MD frames, before returning the averages"""
    DOT = np.zeros(size)
    for t in range(size):
        DOT[t] = np.dot(m1[t], m2[t])
    return np.mean(DOT)

def correlate(residues):
    """Calculates the DCC matrix"""
    sorted_residues = sorted(residues.keys())
    num_trajectories = len(residues[sorted_residues[0]])
    num_residues = len(residues)
    correlation = np.ones((num_residues, num_residues))
    for a, key_a in enumerate(sorted_residues):
        i = residues[key_a]
        resI = np.array(i)
        meanI = np.tile((np.mean(resI, 0)), (num_trajectories, 1))
        idelta = resI - meanI
        magnitudeI = math.sqrt(mean_dot(idelta, idelta, num_trajectories))
        for b, key_b in enumerate(sorted_residues):
            if a < b:
                j = residues[key_b]
                resJ = np.array(j)
                meanJ = np.tile((np.mean(resJ, 0)),(num_trajectories, 1))
                jdelta = resJ - meanJ
                magnitudeJ = math.sqrt(mean_dot(jdelta, jdelta, num_trajectories))
                meanDotIJ = mean_dot(idelta, jdelta, num_trajectories)
                magProd = magnitudeI * magnitudeJ
                corr = meanDotIJ/magProd
                correlation[a, b] = corr
                correlation[b, a] = corr
    return correlation

def plot_map(correlation, title, output_prefix):
    """Plots and saves the heat map"""
    M = np.array(correlation)
    ax = plt.subplots()[1]
    heatmap = ax.pcolor(M, cmap="RdBu_r", vmin=-1, vmax=1)
    fig = plt.gcf()
    ax.set_frame_on(False)
    ax.grid(False)
    plt.xticks(rotation=90)
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(8)
    for t in ax.yaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(8)
    plt.title(title, fontsize=16)
    plt.xlabel('Residue index', fontsize=12)
    plt.ylabel("Residue index", fontsize=12)
    plt.colorbar(heatmap, orientation="vertical")
    plt.savefig('%s.png' % output_prefix, dpi=300)
    plt.close('all')

def write_correlation(correlation, output_prefix):
    """Saves the correlation matrix"""
    with open("%s.txt" % output_prefix, "w") as w:
        rows = correlation.shape[0]
        cols = correlation.shape[1]
        for r in range(rows):
            for c in range(cols):
                w.write('%s ' % str(correlation[r,c]))
            w.write('\n')

def main(args):
    """The main function"""
    log.info("Preparing a trajectory matrix...\n")
    traj_matrix = parse_traj(args.trajectory, args.topology, args.step,
            selected_atoms=args.select_atoms, lazy_load=args.lazy_load)
    log.info("Correlating...\n")
    correlation = correlate(traj_matrix)
    log.info("Plotting heat map...\n")
    plot_map(correlation, args.title, args.prefix)
    write_correlation(correlation, args.prefix)

log = Logger()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Calculates Dynamic Cross Correlation from MD data")
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology",
                        help="Reference PDB file (must contain \
                        the same number of atoms as the trajectory)")
    parser.add_argument("--step",
                        help="Size of the step to take when iterating \
                        the the trajectory frames", type=int)
    parser.add_argument("--select_atoms",
                        help="Comma-separated list of atoms (without spaces) \
                        to use for DCC calculation. E.g. 'CA,P'. The backbone \
                        phosphorus 'P' atom can be a good choice for \
                        representing nucleotides. (default=CA)",
                        type=str, default="CA")
    parser.add_argument("--lazy-load", action='store_true', default=False,
                        help="Iterate through trajectory, loading one \
                        frame into memory at a time (memory-efficient \
                        for large trajectories)")
    parser.add_argument("--title", help="Title for heatmap", default="Protein")
    parser.add_argument("--prefix", default="correlation",
                        help="Prefix for output files")
    log.info("Using Python {}\n".format(sys.version))
    CLI(parser, main, log)
