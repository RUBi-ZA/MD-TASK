#!/usr/bin/env python
#
# Calculate coordination-propensity for all C-Beta/alpha residue pairs 
# 
#
# Script distributed under GNU GPL 3.0
#
# Author: David Penkler
# Date: 09-08-2018

import os, sys, argparse

import numpy 

import mdtraj as md

from math import sqrt

import cp_utils
from lib.cli import CLI
from lib.utils import Logger

def load_trajectory(traj, totalframes, totalres):
    '''Method for loading a trajectory into memory'''
    trajectory = numpy.zeros((totalframes, totalres*3))
    
    for row, frame in enumerate(traj):
        top = frame.topology
        
        col = 0
        for atom_index, atom in enumerate(top.atoms):
            if atom.name == "CA":
                trajectory[row,col:col+3] = frame.xyz[0,atom_index]*10
                col += 3
     
    return trajectory
    
    
def calc_dist(trajectory,totalframes,totalres):
    '''Method for calculating inter-residue distances for every residue pair'''
    distance_matrix = numpy.zeros((totalframes,totalres,totalres))
    print ("___")
    print (str(totalframes))
    print (str(totalres))
    print (str(totalres))
    print ("___")
    for frame in range(0,totalframes):
        for res in range(0,totalres*3,3):
            res_coord = trajectory[frame,res:res+3]
            for other_res in range(0,totalres*3,3):   
                other_coord = trajectory[frame,other_res:other_res+3]
                print ("___")
                print (str(frame))
                print (str(res))
                print (str(other_res))
                print ("___")
                distance_matrix[frame,res/3,other_res/3] = sqrt(((other_coord[0]-res_coord[0])**2)+((other_coord[1]-res_coord[1])**2)+((other_coord[2]-res_coord[2])**2))
    del trajectory
    return distance_matrix


def calculate_cp(distance_matrix,totalframes,totalres):
    '''Method for calculating the coordination propensity for every residue pair'''
    CP = numpy.zeros((totalframes,totalres,totalres))
    for frame in range(0,totalframes):
        for res in range(0,totalres):
            for other_res in range(0,totalres):
                distance = distance_matrix[frame,res,other_res]
                average_distance = numpy.mean(distance_matrix[:,res,other_res])
                cp = (distance-average_distance)**2
                CP[frame,res,other_res] = cp
    CP = numpy.mean(CP,axis=0)
    return CP


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def main(args):
    initial = md.load_frame(args.trajectory, 0, top=args.topology)
    log.info(" -Loading trajectory...")

    if args.num_frames:
        totalframes = args.num_frames
        traj = MDIterator(args.trajectory, top=args.topology, stride=args.step)
    else:
        traj = md.load(args.trajectory, top=args.topology)[::args.step]
        totalframes = traj.n_frames

    totalres = initial.n_residues

    log.info('\tTotal number of frames = %d\n\tNumber of residues = %d' % (totalframes, totalres))

    trajectory = load_trajectory(traj, totalframes, totalres)

    log.info('\tFinal trajectory matrix size: %s' % str(trajectory.shape))
    del traj

    log.info(' -Calculating average distance matrix...')
    distance_matrix = calc_dist(trajectory,totalframes,totalres)
    
    avg_distance = numpy.mean(distance_matrix,axis=0)
    
    log.info(' -Calculating coordination propensity matrix...')
    
    cp_matrix = calculate_cp(distance_matrix,totalframes,totalres)
    numpy.save('%s_cp_matrix.npy' % args.prefix, cp_matrix)
    
    log.info(' -Plotting coordination propensity heatmap...')
    cp_utils.plot_map(cp_matrix,'Coordination propensity', '%s_heatmap' % args.prefix, axes=False)

log = Logger()

if __name__ == "__main__":

    #parse cmd arguments
    parser = argparse.ArgumentParser()
	
    #custom arguments
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Topology PDB file (required if trajectory does not contain topology information)")
    parser.add_argument("--step", help="Size of step when iterating through trajectory frames", default=1, type=int)
    parser.add_argument("--num-frames", help="The number of frames in the trajectory (provides improved performance for large trajectories that cannot be loaded into memory)", type=int, default=None)

    parser.add_argument("--prefix", help="Prefix for output files (default: cp)", default="cp")

    #args = parser.parse_args()
	
    #run script
    #main(args)
	
    CLI(parser, main, log)
