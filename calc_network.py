#!/usr/bin/env python3
# Script distributed under GNU GPL 3.0
# Authors: David Brown and Olivier Sheik Amamuddy
# Date: 23-06-2020

"""Calculate Dynamic Residue interaction Networks from MD trajectories"""

import os
import sys
import argparse
from multiprocessing import Pool as procs

import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

from lib.cli import CLI
from lib.utils import Logger
from lib.trajectory import load_trajectory, calc_distance

def construct_graph(frame, ligands=None, prefix="frame", threshold=6.7,
                    save_graphs=True):
    """Constructs the  network graph for a given MD frame"""
    atom_filter = "(name CB and protein) or (name CA and resname GLY)"
    if ligands:
        ligands = ligands.split(",")
        for ligand in ligands:
            arr = ligand.split(":")
            atom_filter += " or (name {} and resname {})".format(arr[1], arr[0])
    atoms = frame.topology.select(atom_filter)
    nodes_range = len(atoms)
    nodes = range(0, nodes_range)
    edges = []
    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            dist = calc_distance(frame, atoms[i], atoms[j]) * 10
            if dist < threshold:
                edges.append((i, j))
    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)
    if save_graphs:
        nx.write_gml(protein_graph, "{}_graph.gml".format(prefix))
        nx.write_graphml(protein_graph, "{}_graph.graphml".format(prefix))
    return protein_graph

def calc_shortest_path(protein_graph, prefix=None, generate_plots=True):
    """Calculates average shortest path for a given protein graph"""
    num_nodes = len(protein_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)
    path_dict = nx.all_pairs_shortest_path_length(protein_graph)
    path_dict = dict(path_dict)
    dj_path_matrix = np.zeros((num_nodes, num_nodes))
    for i in range(num_nodes):
        for j in range(num_nodes):
            try:
                dj_path_matrix[i, j] = path_dict[i][j]
            except KeyError as _:
                raise nx.exception.NetworkXNoPath(
                    "\nERROR::type=orphan_node:message="
                    "No link between {} and {}:exception=str(_)\n".format(i, j))
    avg_l_per_node = np.sum(dj_path_matrix, axis=0)/(num_nodes - 1)
    avg_l_per_node = pd.Series(avg_l_per_node)
    if generate_plots:
        plt.plot(nodes_axis, avg_l_per_node)
        plt.title("%s L" % prefix, fontsize=18)
        plt.xlabel('Node Indices', fontsize=16)
        plt.ylabel('L', fontsize=16)
        plt.savefig("%s_L.png" % prefix, dpi=300, bbox_inches='tight')
        plt.close()
    avg_l_per_node.to_csv("{}_avg_L.dat".format(prefix), header=False, index=False)
    return avg_l_per_node

def calc_centrality(protein_graph, centrality_function, frameidx=None, 
        label=None, prefix=None, generate_plots=False, centrality_kws=None):
    """Generic function to compute a centrality metric for a frame"""
    if  centrality_kws is None:
        centrality_kws = {}
    try:
        centrality = centrality_function(protein_graph, **centrality_kws)
    except nx.exception.NetworkXNoPath as _:
        log.error("type=orphan_node:frame={}:message=str(_)."
                  "Try increasing the threshold.\n".format(frameidx+1) )
        sys.exit(1)
    centrality = pd.Series(centrality)
    if centrality.sum() == 0:
        log.error("type=unconnected_graph:frame={}:message=str(_)."
                  "Try increasing the threshold.\n".format(frameidx+1) )
        sys.exit(1)

    num_nodes = len(protein_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)
    if generate_plots:
        plt.plot(nodes_axis, centrality)
        plt.title("{} {}".format(prefix, label), fontsize=18)
        plt.xlabel("Node indices", fontsize=16)
        plt.ylabel("{}".format(label.upper()), fontsize=16)
        plt.tight_layout()
        plt.savefig("{}.png".format(prefix), dpi=300, bbox_inches='tight')
        plt.close()
    centrality = pd.Series(centrality)
    centrality.to_csv("{}_{}.dat".format(prefix, label), header=False, index=False)
    return centrality

def calc_centralities(traj, traj_name, total_frames, args):
    """Calculates centralities across all MD frames"""
    log.info("Calculating centralities...\n")
    for current, frame in enumerate(traj):
        try:
            if total_frames:
                log.info("Progress: {}/{}\r".format(current+1, total_frames))
            else:
                log.info("Progress: {} frames completed\r".format(current+1))
            basename = ".".join(traj_name.split(".")[:-1])
            prefix = "{}_{}".format(basename,frame.time[0])
            protein_graph = construct_graph(frame, args.ligands, prefix, args.threshold,
                                            args.save_graphs)
            if args.calc_BC:
                calc_centrality(protein_graph, centrality_function=nx.betweenness_centrality,
                                frameidx=current, prefix=prefix, label="BC",
                                generate_plots=args.generate_plots)
            if args.calc_DC:
                calc_centrality(protein_graph,
                                centrality_function=nx.degree_centrality,
                                frameidx=current, prefix=prefix, label="DC",
                                generate_plots=args.generate_plots)
            if args.calc_CC:
                calc_centrality(protein_graph,
                                centrality_function=nx.closeness_centrality,
                                frameidx=current, prefix=prefix, label="CC",
                                generate_plots=args.generate_plots)
            if args.calc_EC:
                calc_centrality(protein_graph,
                                centrality_function=nx.eigenvector_centrality_numpy,
                                frameidx=current, prefix=prefix, label="EC",
                                generate_plots=args.generate_plots)
            if args.calc_ECC:
                calc_centrality(protein_graph, centrality_function=nx.eccentricity,
                                frameidx=current, prefix=prefix, label="ECC",
                                generate_plots=args.generate_plots)
            if args.calc_PR:
                calc_centrality(protein_graph, centrality_function=nx.pagerank,
                                frameidx=current, prefix=prefix, label="PR",
                                generate_plots=args.generate_plots)
            if args.calc_katz:
                calc_centrality(protein_graph, centrality_function=nx.katz_centrality,
                                frameidx=current, prefix=prefix, label="katz",
                                generate_plots=args.generate_plots)
            if args.calc_L:
                calc_centrality(protein_graph,
                                centrality_function=calc_shortest_path,
                                frameidx=current, prefix=prefix, label="L",
                                generate_plots=args.generate_plots)
        except Exception as _:
            log.error("type=general:frame={}:message={str(_)}\n".format(current + 1))

def main(args):
    """Program start"""
    if not any([args.calc_BC, args.calc_L, args.calc_DC, args.calc_ECC,
                args.calc_EC, args.calc_CC, args.calc_katz, args.calc_PR]):
        log.error("At least one of the network metrics (e.g. --calc-BC) must be set.")
        sys.exit(1)

    results = pd.DataFrame()
    traj_name = os.path.basename(args.trajectory)
    traj, total_frames = load_trajectory(args.trajectory, args.topology,
                                         args.step, args.lazy_load)
    calc_centralities(traj, traj_name, total_frames, args)

log = Logger()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate Dynamic Residue Interaction Networks (DRNs) from \
                         MD trajectories")
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--topology",
                        help="Topology PDB file (required if trajectory does \
                                not contain topology information)")
    parser.add_argument("--ligands", default=None,
                        help="Specify any ligands that should be included in the network")
    parser.add_argument("--threshold", default=6.7, type=float,
                        help="Maximum distance threshold in Angstroms when \
                                constructing graph (default: 6.7)")
    parser.add_argument("--step", default=1, type=int,
                        help="Size of step when iterating through trajectory \
                              frames")
    parser.add_argument("--generate-plots", default=False,
                        help="Generate figures/plots", action='store_true')
    parser.add_argument("--calc-L", action='store_true', default=False,
                        help="Calculate average L")
    parser.add_argument("--calc-BC", action='store_true', default=False,
                        help="Calculate average BC")
    parser.add_argument("--calc-DC", action='store_true', default=False,
                        help="Calculate average degree")
    parser.add_argument("--calc-CC", action='store_true', default=False,
                        help="Calculate average closeness")
    parser.add_argument("--calc-ECC", action='store_true', default=False,
                        help="Calculate average eccentricity")
    parser.add_argument("--calc-EC", action='store_true', default=False,
                        help="Calculate average eigenvector centrality")
    parser.add_argument("--calc-katz", action='store_true', default=False,
                        help="Calculate average Katz centrality")
    parser.add_argument("--calc-PR", action='store_true', default=False,
                        help="Calculate average PageRank")
    parser.add_argument("--save-graphs", action='store_true', default=False,
                        help="Saves network graphs in graphml and gml formats \
                              (default: disabled)")
    parser.add_argument("--lazy-load", action='store_true', default=False,
                        help="Read frames in a memory efficient way - use for \
                              big trajectories)")
    CLI(parser, main, log)
