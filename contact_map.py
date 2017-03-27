#!/usr/bin/env python
#
# Generate weighted contact map around a given residue
#
# Script distributed under GNU GPL 3.0
# 
# Author: Olivier Sheik Amamuddy
# Date: 17-11-2016

import os, sys, argparse

import numpy as np
import subprocess as sp
import mdtraj as md

from datetime import datetime

from lib.utils import *



def write_rscript(script_name, nframes, csv_file, contact_map):
    rscript = """
library('igraph')
setwd("%s")

plot_ = function(grph, edge_weights, layout, title){
    plot(grph,
        edge.width=edge_weights*10,
        edge.label=edge_weights,
        edge.label.color="black",
        edge.label.cex=0.5,
        vertex.color="darkolivegreen3",
        vertex.size=30,
        vertex.frame.color=NA,
        vertex.label.cex=.65,
        vertex.label.color="black",
        vertex.label.family="serif",
        main=title, layout=layout
    )
}

draw.network = function(dat){
    dat.fmt = cbind(gsub("_",".", dat[,1]), gsub("_",".", dat[,2]))
    g = graph_from_edgelist(dat.fmt, directed = 'FALSE')
    adjMat = as_adjacency_matrix(g)
    adjMat = (adjMat/%d)
    g = graph_from_adjacency_matrix(adjMat, weighted =T, mode='undirected')
    edge_weights = round(E(g)$weight, 3)

    title = gsub("_network.csv","", "%s")
    title = paste0(title, '_contact network')

    l = layout_with_dh(g, maxiter=50)
    plot_(grph=g, edge_weights=edge_weights, layout=l, title=title)
}

pdf("%s", paper='a4r', width=12, height=12)
par(mfrow=c(1,1), cex.main=.8)

dat = read.csv("%s", header=F)
draw.network(dat=dat)

#Focusing on the variants
dev.off()
""" % (os.getcwd(), nframes, csv_file, contact_map, csv_file)

    with open(script_name, "w") as f:
        f.write(rscript)


def main(args):
    traj_path = args.trajectory
    topology = args.topology
    residue  = args.residue.upper()
    cutoff = args.threshold / 10
    prefix = residue

    chainChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    edges = []
    
    log("Loading trajectory...")
    
    traj = md.load(traj_path, top=topology)[::args.step]

    # Get CA and CB atom indices
    atom_indices = [atom.index for atom in traj.topology.atoms if ((atom.name == "CB" and atom.residue.name != "GLY") or \
        (atom.name == "CA" and atom.residue.name == "GLY"))]

    residues = list(map(lambda x: str(x), traj.top.residues))
    nframes = len(traj)

    log("Calculating network around %s..." % residue)
    
    for frame in traj:

        # Takes a frame and updates a list of contact
        for aaidx, atom1_idx in enumerate(atom_indices):

            if residues[aaidx] == residue:
                atom1          = list(frame.top.atoms)[atom1_idx]
                atom1_chain    = chainChar[atom1.residue.chain.index]
                atom1_resname  = atom1.residue.name
                atom1_resid    = atom1.residue.resSeq

                for atom2_idx in atom_indices:
                    if atom1_idx != atom2_idx:
                        # Calculate distances and create edge
                        distance = np.linalg.norm(frame.xyz[0,atom1_idx] - frame.xyz[0,atom2_idx])

                        if distance < cutoff:
                            atom2         = list(frame.top.atoms)[atom2_idx]
                            atom2_chain   = chainChar[atom2.residue.chain.index]
                            atom2_resname = atom2.residue.name
                            atom2_resid   = atom2.residue.resSeq
    
                            # Write contact
                            contact = "{}{}_{},{}{}_{}".format(
                            atom1_resname, atom1_resid, atom1_chain,
                            atom2_resname, atom2_resid, atom2_chain)
                            edges.append(contact)

                break #Stop looking for other positions

    #Write one variant at a time
    csv_file = "%s_network.csv" % prefix
    contact_map = "%s_contact_map.pdf" % prefix

    log("Writing network to %s..." % csv_file)

    with open(csv_file, 'w') as f_handle:
        edges = "\n".join(edges)
        f_handle.write(edges)

    log("Generating contact map: %s..." % contact_map)
    
    script_name = "rscript.R"
    write_rscript(script_name, nframes, csv_file, contact_map)
        
    sp.call("R CMD BATCH %s" % script_name, shell=True)

    if os.path.exists(contact_map):
        os.remove(script_name)
    else:
        r_out = "rscript.out"
        log("Contact map not generated. See contents of %s for details..." % r_out)



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
    parser.add_argument("--residue", help="The residue that the contact map will be built around (e.g. THR405)")
    parser.add_argument("--threshold", help="Maximum distance threshold in Angstroms when constructing graph (default: 6.7)", default=6.7, type=float)
    parser.add_argument("--prefix", help="Prefix for output file", default="residue_contacts")
    parser.add_argument("--step", help="Size of step when iterating through trajectory frames", default=1, type=int)
    
    args = parser.parse_args()
    
    #set up logging
    silent = args.silent
    
    if args.log_file:
        stream = open(args.log_file, 'w')
    
    start = datetime.now()
    log("Started at: %s\n" % str(start))

    #run script
    main(args)

    end = datetime.now()
    time_taken = format_seconds((end - start).seconds)
    
    log("\nCompleted at: %s\n" % str(end))
    log("- Total time: %s\n" % str(time_taken))
    
    #close logging stream
    stream.close()
    
 
