#!/usr/bin/env python
#
# Perform PRS calculations given and MD trajectory and a final state
# co-ordinate file
#
# Script distributed under GNU GPL 3.0
#
# Author: David Penkler
# Date: 17-11-2016

import sys
import argparse
from math import log10, floor, sqrt
import numpy
import mdtraj as md
from lib import sdrms
from lib.cli import CLI
from lib.utils import Logger
from lib.trajectory import load_trajectory

def round_sig(num, sig=2):
    """Rounds to x significant numbers"""
    return round(num, sig-int(floor(log10(num)))-1)

def coarsegrain_topology(topology_filename, sele="name CA", save_xyz=False):
    """Coarsegrain topology
    Returns:
    outfilename, coarsegrained mdtraj topology in Ang, if saved
    """
    topology = md.load(topology_filename)
    selection = topology.top.select(sele)
    topology.atom_slice(selection, inplace=True)
    if save_xyz:
        outfilename = "{}.xyz".format(topology_filename[:-4])
        topology.save_xyz(outfilename)
        log.info("Wrote initial co-ordinate file {}\n".format(outfilename))
        return outfilename, topology
    return topology

def trajectory_to_array(traj, totalframes, totalres, selection="name CA"):
    """Reshape trajectory to nframes x (3*n_res) and return trajectory in Ang"""
    selection = traj.top.select(selection)
    traj.atom_slice(selection, inplace=True)
    trajectory = traj.xyz.reshape(totalframes, totalres*3)
    return trajectory*10

def align_frame(reference_frame, alternative_frame, aln=False, mask=None):
    """"Align each frame of the trajectory to a reference"""
    totalres = reference_frame.shape[0]
    if aln:
        return sdrms.superpose3D(alternative_frame.reshape(totalres, 3),
                                 reference_frame, refmask=mask,
                                 targetmask=mask)[0].reshape(1, totalres*3)[0]
    return sdrms.superpose3D(alternative_frame.reshape(totalres, 3),
                             reference_frame)[0].reshape(1, totalres*3)[0]

def calc_rmsd(reference_frame, alternative_frame, aln=False, mask=None):
    """Calculate CA atom RMSD"""
    if aln:
        return sdrms.superpose3D(alternative_frame,
                                 reference_frame, refmask=mask,
                                 targetmask=mask)[1]
    return sdrms.superpose3D(alternative_frame, reference_frame)[1]

def main(args):
    if not args.final:
        log.error("a final co-ordinate file must be supplied via the --final argument\n")
        sys.exit(1)
    # Use topology as initial conformation
    args.initial, initial = coarsegrain_topology(args.topology, save_xyz=True)
    args.final, final = coarsegrain_topology(args.final, save_xyz=True)
    # Disabling aln due to missing mask
    args.aln = False
    mask = None
    # Convert coorsd to Angstrom
    initial.xyz *= 10
    final.xyz *= 10
    log.info("Loading trajectory...\n")
    if args.num_frames:
        traj, totalframes = load_trajectory(args.trajectory, args.topology,
                                            args.step, True)
        totalframes = args.num_frames
    else:
        traj, totalframes = load_trajectory(args.trajectory, args.topology,
                                            args.step, False)
    totalres = initial.n_residues
    log.info('- Total number of frames = %d\n- Number of residues = %d\n' % (totalframes, totalres))
    trajectory = trajectory_to_array(traj, totalframes, totalres)
    log.info('- Final trajectory matrix size: %s\n' % str(trajectory.shape))
    del traj
    log.info("Aligning trajectory frames...\n")
    aligned_mat = numpy.zeros((totalframes,3*totalres))
    frame_0 = trajectory[0].reshape(totalres, 3)

    for frame in range(0, totalframes):
        aligned_mat[frame] = align_frame(frame_0, trajectory[frame], args.aln)
    del trajectory
    log.info("- Calculating average structure...\n")
    average_structure_1 = numpy.mean(aligned_mat, axis=0).reshape(totalres, 3)
    log.info("- Aligning to average structure...\n")
    for i in range(0, 10):
        for frame in range(0, totalframes):
            aligned_mat[frame] = align_frame(average_structure_1,
                                             aligned_mat[frame],
                                             args.aln)
        average_structure_2 = numpy.average(aligned_mat, axis=0).reshape(totalres, 3)
        rmsd = calc_rmsd(average_structure_1, average_structure_2, args.aln)
        log.info('   - %s Angstroms from previous structure\n' % str(rmsd))
        average_structure_1 = average_structure_2
        del average_structure_2
        if rmsd <= 0.000001:
            for frame in range(0, totalframes):
                aligned_mat[frame] = align_frame(average_structure_1,
                                                 aligned_mat[frame],
                                                 args.aln)
            break
    log.info("Calculating difference between frame atoms and average atoms...\n")
    meanstructure = average_structure_1.reshape(totalres*3)
    del average_structure_1
    log.info('- Calculating R_mat\n')
    R_mat = numpy.zeros((totalframes, totalres*3))
    for frame in range(0, totalframes):
        R_mat[frame,:] = (aligned_mat[frame,:]) - meanstructure
    log.info('- Transposing\n')
    RT_mat = numpy.transpose(R_mat)
    RT_mat = numpy.mat(RT_mat)
    R_mat = numpy.mat(R_mat)
    log.info('- Calculating corr_mat\n')
    corr_mat = (RT_mat * R_mat)/ (totalframes-1)
    numpy.savetxt("corr_mat.txt", corr_mat)
    del aligned_mat
    del meanstructure
    del R_mat
    del RT_mat
    log.info('Reading initial and final PDB co-ordinates...\n')
    initial = initial.xyz[0]
    final = final.xyz[0]
    if initial.shape[0] != final.shape[0]:
        log.error("Unequal number of CA atoms in initial and final structures. Check provided structures.")
        sys.exit(1)
    log.info('Calculating experimental difference between initial and final co-ordinates...\n')
    if args.aln:
        log.info("- Using NTD alignment restrictions\n")
        final_alg = sdrms.superpose3D(final, initial, refmask=mask,
                                      targetmask=mask)[0]
    else:
        final_alg = sdrms.superpose3D(final, initial)[0]
    diffE = (final_alg-initial).reshape(totalres*3, 1)
    del final
    del final_alg
    log.info('Implementing perturbations sequentially...\n')
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
                diffP[i,:,s] = ((sdrms.superpose3D(diffP[i,:,s].reshape(totalres, 3),
                                                   initial, refmask=mask,
                                                   targetmask=mask)[0].reshape(1, totalres*3))[0]) - initial_trans[0]
            else:
                diffP[i,:,s] = ((sdrms.superpose3D(diffP[i,:,s].reshape(totalres, 3),
                                                   initial)[0].reshape(1, totalres*3))[0]) - initial_trans[0]
            del delF
    del initial_trans
    del initial
    del corr_mat
    log.info("Calculating Pearson's correlations coefficient...\n")
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


log = Logger()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology", required=True,
            help="Topology PDB file for the initial state co-ordinate file (required)")
    parser.add_argument("--step", default=1, type=int, required=True,
            help="Size of step when iterating through trajectory frames")
    parser.add_argument("--final", help="Final state co-ordinate file (required)")
    parser.add_argument("--perturbations", type=int, default=250,
            help="Number of perturbations (default: 250)")
    parser.add_argument("--num-frames", type=int, default=None,
            help="The number of frames in the trajectory (provides improved performance for large trajectories that cannot be loaded into memory)")
    parser.add_argument("--aln", help="Restrict N-Terminal alignment (not implemented)", action="store_true")
    parser.add_argument("--prefix", default="result",
            help="Prefix for CSV output file (default: result)")
    CLI(parser, main, log)
