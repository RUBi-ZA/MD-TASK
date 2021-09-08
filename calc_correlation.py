#!/usr/bin/env python
#
# Calculate correlation in MD trajectory
#
# Script distributed under GNU GPL 3.0
#
# Author: Caroline Ross
# Date: 17-11-2016

from matplotlib import cm

import numpy as np

from lib.cli import CLI
from lib.utils import Logger
from lib.trajectory import load_trajectory

import argparse, math, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_traj(traj, topology=None, step=1, selected_atoms=["CA"], lazy_load=False):
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
                    residues[res] =  [co_ords]

    return residues



def mean_dot(m1, m2, size):
    DOT = np.zeros(size)

    for t in range(size):
        DOT[t] = np.dot(m1[t],m2[t])

    return np.mean(DOT)



def correlate(residues):

    sorted_residues = sorted(residues.keys())

    num_trajectories = len(residues[sorted_residues[0]])
    num_residues = len(residues)

    correlation = np.zeros((num_residues, num_residues))

    for a, key_a in enumerate(sorted_residues):
        i = residues[key_a]
        resI = np.array(i)
        meanI = np.tile((np.mean(resI, 0)),(num_trajectories, 1))
        idelta = resI - meanI;
        magnitudeI = math.sqrt(mean_dot(idelta, idelta, num_trajectories))

        for b, key_b in enumerate(sorted_residues):
            j = residues[key_b]
            resJ = np.array(j)
            meanJ = np.tile((np.mean(resJ, 0)),(num_trajectories, 1))
            jdelta = resJ - meanJ
            magnitudeJ = math.sqrt(mean_dot(jdelta, jdelta, num_trajectories))

            meanDotIJ = mean_dot(idelta, jdelta, num_trajectories)
            magProd = magnitudeI * magnitudeJ

            correlation[a,b] = meanDotIJ/magProd

    return correlation


def plot_map(correlation, title, output_prefix):
    M = np.array(correlation)

    ax = plt.subplots()[1]
    colors = [('white')] + [(cm.jet(i)) for i in range(40,250)]

    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(M, cmap=new_map, vmin=-1, vmax=1)

    fig = plt.gcf()
    ax.set_frame_on(False)
    ax.grid(False)

    plt.xticks(rotation=90)

    # Turn off all the ticks
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
    plt.xlabel('Residue Index', fontsize=12)
    plt.ylabel("Residue Index", fontsize=12)

    plt.colorbar(heatmap, orientation="vertical")
    plt.savefig('%s.png' % output_prefix, dpi=300)
    plt.close('all')



def print_correlation(correlation, output_prefix):
    with open("%s.txt" % output_prefix, "w") as w:
        rows = correlation.shape[0]
        cols = correlation.shape[1]

        for r in range(rows):
            for c in range(cols):
                w.write('%s ' % str(correlation[r,c]))
            w.write('\n')



def main(args):
    log.info("Preparing a trajectory matrix...\n")
    traj_matrix = parse_traj(args.trajectory, args.topology, args.step, lazy_load=args.lazy_load)

    log.info("Correlating...\n")
    correlation = correlate(traj_matrix)

    log.info("Plotting heat map...\n")
    plot_map(correlation, args.title, args.prefix)
    print_correlation(correlation, args.prefix)


log = Logger()

if __name__ == "__main__":

    #parse cmd arguments
    parser = argparse.ArgumentParser()

    #custom arguments
    parser.add_argument("--trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Referencce PDB file (must contain the same number of atoms as the trajectory)")
    parser.add_argument("--step", help="Size of the step to take when iterating the the trajectory frames", type=int)
    parser.add_argument("--lazy-load", help="Iterate through trajectory, loading one frame into memory at a time (memory-efficient for large trajectories)", action='store_true', default=False)

    parser.add_argument("--title", help="Title for heatmap", default="Protein")
    parser.add_argument("--prefix", help="Prefix for output files", default="correlation")

    CLI(parser, main, log)
