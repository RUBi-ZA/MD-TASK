#!/usr/bin/env python
#
# Calculate residue interaction network for frames in an MD trajectory and 
# determine betweenness centrality and average shortest path for residues in
# the frames
#
# Script distributed under GNU GPL 3.0
# 
# Author: David Brown
# Date: 17-11-2016


#from prody import *
import mdtraj as md

#from natsort import natsorted

from lib.utils import *

from datetime import datetime

import numpy as np
import networkx as nx

import os, sys, argparse, math, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt



def calcDistance(frame, index1, index2):
    atom1 = frame.xyz[0, index1]
    atom2 = frame.xyz[0, index2]
    
    dist = math.sqrt((atom2[0] - atom1[0])**2 + (atom2[1] - atom1[1])**2 + (atom2[2] - atom1[2])**2)
    
    return abs(dist)
    


def construct_graph(frame, ligands=None, prefix="frame", threshold=6.7, save_graph=True):
    filter = "(name CB and protein) or (name CA and resname GLY)"
    if ligands:
        ligands = ligands.split(",")
        
        for ligand in ligands:
            arr = ligand.split(":")
            filter += " or (name %s and resname %s)" % (arr[1], arr[0])
    
    atoms = frame.topology.select(filter)
    
    nodes_range = len(atoms)
    
    nodes = range(0, len(atoms))
    edges = []
    
    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            dist = calcDistance(frame, atoms[i], atoms[j]) * 10
            if dist < threshold:
                edges.append((i, j))

    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)
    
    if save_graph:
        nx.write_gml(protein_graph, "%s_graph.gml" % prefix)
        nx.write_graphml(protein_graph, "%s_graph.graphml" % prefix)
    
    return protein_graph


def calc_shortest_path(protein_graph, prefix, generate_plots=True):
    num_nodes = len(protein_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)
    
    dj_path_matrix = np.zeros((num_nodes, num_nodes))
    
    for i in range(num_nodes):
        for j in range(num_nodes):
            dj_path_matrix[i,j] = nx.dijkstra_path_length(protein_graph, i, j)
    
    np.savetxt("%s_L.dat" % prefix, dj_path_matrix)
            
    avg_L_per_node = np.sum(dj_path_matrix, axis=0)/(num_nodes - 1)
    
    if generate_plots:
        plt.plot(nodes_axis, avg_L_per_node)
        plt.title("%s L" % prefix, fontsize=18)
        plt.xlabel('Node Indices', fontsize=16)
        plt.ylabel('L', fontsize=16)
        plt.savefig("%s_L.png" % prefix, dpi=300, bbox_inches='tight')
        plt.close()
    
    avg_L_per_node = avg_L_per_node.reshape(1, num_nodes)
    np.savetxt("%s_avg_L.dat" % prefix, avg_L_per_node)
    
    return dj_path_matrix


def calc_BC(protein_graph, prefix, generate_plots=True):
    bc = nx.betweenness_centrality(protein_graph, normalized=True)
    bc = np.asarray(list(bc.values()))
    
    num_nodes = len(protein_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)
    
    if generate_plots:
        plt.plot(nodes_axis, bc)
        plt.title("%s BC" % prefix, fontsize=18)
        plt.xlabel('Node Indices', fontsize=16)
        plt.ylabel('BC', fontsize=16)
        plt.savefig("%s_BC.png" % prefix, dpi=300, bbox_inches='tight')
        plt.close()
    
    bc = bc.reshape(1, num_nodes)
    np.savetxt("%s_bc.dat" % prefix, bc)
    
    return bc


def main(args):
    global stream
    
    traj = md.load(args.trajectory, top=args.topology)[::args.step]
    
    total_frames = len(traj)
    
    traj_name = os.path.basename(args.trajectory)
    
    for current, frame in enumerate(traj):
        log("Progress: %d/%d\r" % (current + 1, total_frames)) 
        
        prefix = "%s_%d" % (".".join(traj_name.split(".")[:-1]), frame.time)
        
        pg = construct_graph(frame, args.ligands, prefix, args.threshold, args.discard_graphs)
        if args.calc_L:
            calc_shortest_path(pg, prefix, args.generate_plots)
        if args.calc_BC:
            calc_BC(pg, prefix, args.generate_plots)
        
        stream.flush()


silent = False
stream = sys.stdout

def log(message):
    global silent
    global stream
    
    if not silent:
        stream.write(str(message))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    #standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    
	#custom arguments 
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology", help="Topology PDB file (required if trajectory does not contain topology information)")
    parser.add_argument("--ligands", help="Specify any ligands that should be included in the network", default=None)
    parser.add_argument("--threshold", help="Maximum distance threshold in Angstroms when constructing graph (default: 6.7)", default=6.7, type=float)
    
    parser.add_argument("--step", help="Size of step when iterating through trajectory frames", default=1, type=int)
    parser.add_argument("--generate-plots", help="Generate figures/plots", action='store_true', default=False)
    parser.add_argument("--calc-L", help="Calculate delta L", action='store_true', default=False)
    parser.add_argument("--calc-BC", help="Calculate delta BC", action='store_true', default=False)
    parser.add_argument("--discard-graphs", help="Discard calculated networks when complete (default: save networks in graphml and gml formats)", action='store_false', default=True)
    
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
    
