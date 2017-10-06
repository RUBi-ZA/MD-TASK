#!/usr/bin/env python
#
# Perform PRS calculations given and MD trajectory and a final state 
# co-ordinate file
#
# Script distributed under GNU GPL 3.0
# 
# Author: David Penkler
# Date: 17-11-2016

import os, sys, argparse 

import numpy, scipy

import mdtraj as md

from scipy.stats.stats import pearsonr
from scipy import linalg as LA

from math import log10, floor, sqrt
from datetime import datetime

from lib import sdrms
from lib.utils import *
 


def round_sig(x, sig=2):
    return round(x,sig-int(floor(log10(x)))-1)


def load_trajectory(traj, totalframes, totalres):
    trajectory = numpy.zeros((totalframes, totalres*3))
    
    for row, frame in enumerate(traj):
        top = frame.topology
        
        col = 0
        for atom_index, atom in enumerate(top.atoms):
            if atom.name == "CA":
                trajectory[row,col:col+3] = frame.xyz[0,atom_index]*10
                col += 3
     
    return trajectory


def align_frame(reference_frame, alternative_frame, aln=False):
    totalres = reference_frame.shape[0]

    if aln:
        return sdrms.superpose3D(alternative_frame.reshape(totalres, 3), reference_frame, refmask=mask, targetmask=mask)[0].reshape(1, totalres*3)[0]
    else:
        return sdrms.superpose3D(alternative_frame.reshape(totalres, 3), reference_frame)[0].reshape(1, totalres*3)[0]


def calc_rmsd(reference_frame, alternative_frame, aln=False):
    if aln:
        return sdrms.superpose3D(alternative_frame, reference_frame, refmask=mask, targetmask=mask)[1]
    else:
        return sdrms.superpose3D(alternative_frame, reference_frame)[1]



def main(args):
    if not args.final:
        log("ERROR: a final co-ordinate file must be supplied via the --final argument")
        sys.exit(1) 
    
    initial = md.load_frame(args.trajectory, 0, top=args.topology)
    if not args.initial:
        args.initial = "initial.xyz"       
        
        log("Generating initial co-ordinate file: %s" % args.initial)
        initial[0].save(args.initial)
    
    
    log("Loading trajectory...")
    
    if args.num_frames:
        totalframes = args.num_frames
        traj = MDIterator(args.trajectory, top=args.topology, stride=args.step)   
    else:
        traj = md.load(args.trajectory, top=args.topology)[::args.step]  
        totalframes = traj.n_frames
    
    totalres = initial.n_residues
    
    log('- Total number of frames = %d\n- Number of residues = %d' % (totalframes, totalres))
    
    trajectory = load_trajectory(traj, totalframes, totalres) 
    
    log('- Final trajectory matrix size: %s' % str(trajectory.shape))    
    del traj
    
    
    log("Aligning trajectory frames...")

    aligned_mat = numpy.zeros((totalframes,3*totalres))
    frame_0 = trajectory[0].reshape(totalres, 3)
    
    for frame in range(0, totalframes):
        aligned_mat[frame] = align_frame(frame_0, trajectory[frame], args.aln)
    
    del trajectory
    
    
    log("- Calculating average structure...")
                    
    average_structure_1 = numpy.mean(aligned_mat, axis=0).reshape(totalres, 3)
        
    
    log("- Aligning to average structure...")
    
    for i in range(0, 10):
        for frame in range(0, totalframes):
            aligned_mat[frame] = align_frame(average_structure_1, aligned_mat[frame], args.aln)
        
        average_structure_2 = numpy.average(aligned_mat, axis=0).reshape(totalres, 3)
        
        rmsd = calc_rmsd(average_structure_1, average_structure_2, args.aln)
        
        log('   - %s Angstroms from previous structure' % str(rmsd))
        
        average_structure_1 = average_structure_2
        del average_structure_2
        
        if rmsd <= 0.000001:
            for frame in range(0, totalframes):
                aligned_mat[frame] = align_frame(average_structure_1, aligned_mat[frame], args.aln)
            break
    
    
    log("Calculating difference between frame atoms and average atoms...")
    
    meanstructure = average_structure_1.reshape(totalres*3)
        
    del average_structure_1
    
    log('- Calculating R_mat')
    R_mat = numpy.zeros((totalframes, totalres*3))
    for frame in range(0, totalframes):
        R_mat[frame,:] = (aligned_mat[frame,:]) - meanstructure  
    
    log('- Transposing')
    
    RT_mat = numpy.transpose(R_mat)
        
    RT_mat = numpy.mat(RT_mat)
    R_mat = numpy.mat(R_mat)
    
    log('- Calculating corr_mat')
    
    corr_mat = (RT_mat * R_mat)/ (totalframes-1)
    numpy.savetxt("corr_mat.txt", corr_mat)  

    del aligned_mat
    del meanstructure
    del R_mat
    del RT_mat
    
    
    log('Reading initial and final PDB co-ordinates...')
        
    initial = numpy.zeros((totalres, 3))
    final = numpy.zeros((totalres, 3))
    
    with open(args.initial, 'r') as initial_lines:
        with open(args.final, 'r') as final_lines:
            
            res_index = 0
            for line_index, initial_line in enumerate(initial_lines):
                final_line = final_lines.readline()
                
                if line_index >= 2 and res_index < totalres:
                    initial_res = initial_line.strip().split()
                    
                    if initial_res[0] == "CA":
                        final_res = final_line.strip().split()
                        
                        initial[res_index,] = initial_res[1:]
                        final[res_index,] = final_res[1:]
                        res_index += 1
           
        
    log('Calculating experimental difference between initial and final co-ordinates...')
    
    if args.aln:
        log("- Using NTD alignment restrictions")
        final_alg = sdrms.superpose3D(final, initial, refmask=mask, targetmask=mask)[0]
    else:
        final_alg = sdrms.superpose3D(final, initial)[0]

    diffE = (final_alg-initial).reshape(totalres*3, 1)

    del final
    del final_alg
    
    
    log('Implementing perturbations sequentially...')
    
    perturbations = int(args.perturbations)
    diffP = numpy.zeros((totalres, totalres*3, perturbations))
    initial_trans = initial.reshape(1, totalres*3)
     
    for s in range(0, perturbations):
        for i in range(0, totalres):
            delF = numpy.zeros((totalres*3))
            f = 2 * numpy.random.random((3, 1)) - 1
            j = (i + 1) * 3
            
            delF[j-3] = round_sig(abs(f[0,0]), 5)* -1 if f[0,0]< 0 else round_sig(abs(f[0,0]), 5)
            delF[j-2] = round_sig(abs(f[1,0]), 5)* -1 if f[1,0]< 0 else round_sig(abs(f[1,0]), 5)
            delF[j-1] = round_sig(abs(f[2,0]), 5)* -1 if f[2,0]< 0 else round_sig(abs(f[2,0]), 5)
            
            diffP[i,:,s] = numpy.dot((delF), (corr_mat))
            diffP[i,:,s] = diffP[i,:,s] + initial_trans[0]
            
            if args.aln:
                diffP[i,:,s] = ((sdrms.superpose3D(diffP[i,:,s].reshape(totalres, 3), initial, refmask=mask, targetmask=mask)[0].reshape(1, totalres*3))[0]) - initial_trans[0]
            else:
                diffP[i,:,s] = ((sdrms.superpose3D(diffP[i,:,s].reshape(totalres, 3), initial)[0].reshape(1, totalres*3))[0]) - initial_trans[0]
            del delF
    
    del initial_trans
    del initial
    del corr_mat
    
    
    log("Calculating Pearson's correlations coefficient...")

    DTarget = numpy.zeros(totalres)
    DIFF = numpy.zeros((totalres, totalres, perturbations))
    RHO = numpy.zeros((totalres, perturbations))

    for i in range(0, totalres):
        DTarget[i] = sqrt(diffE[3*(i+1)-3]**2 + diffE[3*(i+1)-2]**2 + diffE[3*(i+1)-1]**2)
    
    for j in range(0, perturbations):
        for i in range(0, totalres):
            for k in range(0, totalres):
                DIFF[k,i,j] = sqrt((diffP[i, 3*(k+1)-3, j]**2) + (diffP[i, 3*(k+1)-2, j]**2) + (diffP[i, 3*(k+1)-1, j]**2))
    
    del diffP

    for i in range(0, perturbations):
        for j in range(0, totalres):
            RHO[j,i] = numpy.corrcoef(numpy.transpose(DIFF[:,j,i]), DTarget)[0,1]
    
    del DIFF
    del DTarget            
          
    maxRHO = numpy.zeros(totalres)
    for i in range(0, totalres):
        maxRHO[i] = numpy.amax(abs(RHO[i,:]))
        
    numpy.savetxt("%s.csv" % args.prefix, maxRHO, delimiter=",", header=args.prefix)
    
    del maxRHO



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
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Topology PDB file (required if trajectory does not contain topology information)")
    parser.add_argument("--step", help="Size of step when iterating through trajectory frames", default=1, type=int)
    parser.add_argument("--initial", help="Initial state co-ordinate file (default: generated from first frame of trajectory)", default=None)
    parser.add_argument("--final", help="Final state co-ordinate file (must be provided)")
    parser.add_argument("--perturbations", help="Number of perturbations (default: 250)", type=int, default=250)
    parser.add_argument("--num-frames", help="The number of frames in the trajectory (provides improved performance for large trajectories that cannot be loaded into memory)", type=int, default=None)
    parser.add_argument("--aln", help="Restrict N-Terminal alignment", action="store_true")
    
    parser.add_argument("--prefix", help="Prefix for CSV output file (default: result)", default="result")
    
    args = parser.parse_args()
    
    #set up logging
    silent = args.silent
    
    if args.log_file:
        stream = open(args.log_file, 'w')
    
    start = datetime.now()
    log("Started at: %s" % str(start))
    
    #run script
    main(args)
    
    end = datetime.now()
    time_taken = format_seconds((end - start).seconds)
    
    log("Completed at: %s" % str(end))
    log("- Total time: %s" % str(time_taken))
    
    #close logging stream
    stream.close()

        

