#!/usr/bin/env python
#
# Calculate correlation in MD trajectory 
#
# Script distributed under GNU GPL 3.0
# 
# Author: Caroline Ross
# Date: 17-11-2016


import os, sys, argparse, traceback

import math, matplotlib
import matplotlib.pyplot as plt

from matplotlib import cm

import mdtraj as md

import numpy as np



def parse_traj(traj, topology=None, step=1, selected_atoms=["CA"]):    
    traj = md.load(traj, top=topology)    
    
    residues = {}
    
    for frame in traj[::step]:        
        for atom in frame.topology.atoms:
            if atom.name in selected_atoms:
                res = atom.residue.resSeq
                
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

    sorted_residues = sorted(residues.iterkeys())
    
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



def plot_map(correlation, output_prefix):   
    M = np.array(correlation)
    
    fig, ax = plt.subplots()
    colors = [('white')] + [(cm.jet(i)) for i in xrange(40,250)]
    
    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(M, cmap=new_map, vmin=-1, vmax=1)
    
    fig = plt.gcf()
    fig.set_size_inches(20, 20)
    ax.set_frame_on(False)
    ax.grid(False)
    
    plt.xticks(rotation=90)
    
    # Turn off all the ticks
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    
    cbar = plt.colorbar(heatmap, orientation="vertical")
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
    log("Preparing a trajectory matrix...")
    traj_matrix = parse_traj(args.trajectory, args.topology, args.step)
    
    log("Correlating...")
    correlation = correlate(traj_matrix)
    
    log("Plotting heat map...")
    plot_map(correlation, args.prefix)
    print_correlation(correlation, args.prefix)
    
    log("Completed successfully!")



silent = False
stream = sys.stdout

def log(message):
    global silent
    global stream
    
    if not silent:
        print >> stream, message
    


if __name__ == "__main__":
    
    #parse cmd arguments
    parser = argparse.ArgumentParser()
    
    #standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    
    #custom arguments
    parser.add_argument("--trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Referencce PDB file (must contain the same number of atoms as the trajectory)")
    parser.add_argument("--step", help="Size of the step to take when iterating the the trajectory frames", type=int)

    parser.add_argument("--prefix", help="Prefix for output files")
    
    args = parser.parse_args()
    
    #set up logging
    silent = args.silent
    
    if args.log_file:
        stream = open(args.log_file, 'w')
    
    #run script
    main(args)
    
    #close logging stream
    stream.close()




 
