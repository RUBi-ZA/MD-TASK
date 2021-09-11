#!/usr/bin/env python
"""
Calculates various DRN metrics
"""


__author__ = ["Olivier Sheik Amamuddy", "David Brown"]
__copyright__ = "Copyright 2020, Research Unit in Bioinformatics"
__license__ = "GNU GPL 3.0"
__version__ = "1.1"
__email__ = "oliserand@gmail.com"
__status__ = "Production"
__date__ = "7th Oct 2020"

import os
import sys
import argparse
from multiprocessing import Pool

import mdtraj as md
import pandas as pd
import numpy as np
import networkx as nx
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
import matplotlib.pyplot as plt
import matplotlib

from lib.cli import CLI
from lib.utils import Logger
from lib.trajectory import load_trajectory, calc_distance
matplotlib.use('Agg')


def get_selection(atom_filter="(name CB) or (name CA and resname GLY)",
                  ligands=False):
    """Prepares the selection string for use in atom selection from a mdtraj
    frame(s)"""
    if ligands:
        ligands = ligands.split(",")
        for ligand in ligands:
            arr = ligand.split(":")
            atom_filter += " or (name {} and resname {})".format(arr[1], arr[0])
    return atom_filter

def construct_graph(frame, atom_filter=None, prefix="frame", threshold=6.7,
                    save_graphs=True):
    """Constructs the  network graph for a given MD frame"""
    atoms = frame.topology.select(atom_filter)
    nodes_range = len(atoms)
    nodes = range(0, nodes_range)
    edges = []
    threshold /= 10
    for i in range(nodes_range - 1):
        for j in range(i + 1, nodes_range):
            dist = calc_distance(frame, atoms[i], atoms[j])
            if dist < threshold:
                edges.append((i, j))
    p_graph = nx.Graph()
    p_graph.add_nodes_from(nodes)
    p_graph.add_edges_from(edges)
    if save_graphs:
        nx.write_gml(p_graph, "{}_graph.gml".format(prefix))
        nx.write_graphml(p_graph, "{}_graph.graphml".format(prefix))
    return p_graph

def calc_shortest_path(p_graph, prefix=None, generate_plots=True):
    """Calculates average shortest path for a given protein graph"""
    num_nodes = len(p_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)
    path_dict = nx.all_pairs_shortest_path_length(p_graph)
    path_dict = dict(path_dict)
    dj_path_matrix = np.zeros((num_nodes, num_nodes))
    for i in range(num_nodes):
        for j in range(num_nodes):
            try:
                dj_path_matrix[i, j] = path_dict[i][j]
            except KeyError:
                LOG.error("\nERROR::type=orphan_node:message=No link between {0} and {1}:exception={2}\n".format(i, j, nx.exception.NetworkXNoPath))
                raise nx.exception.NetworkXNoPath
    avg_l_per_node = np.sum(dj_path_matrix, axis=0)/(num_nodes - 1)
    np.savetxt("Dij.csv", dj_path_matrix, delimiter=",") #added
    avg_l_per_node = pd.Series(avg_l_per_node)
    if generate_plots:
        plt.plot(nodes_axis, avg_l_per_node)
        plt.title("%s L" % prefix, fontsize=18)
        plt.xlabel('Node Indices', fontsize=16)
        plt.ylabel('L', fontsize=16)
        plt.savefig("%s_avg_L.png" % prefix, dpi=300, bbox_inches='tight')
        plt.close()
    return avg_l_per_node

def calc_centrality(p_graph, c_func, frameidx=None,
                    label=None, prefix=None, generate_plots=False,
                    centrality_kws=None):
    """Generic function to compute a centrality metric for a frame"""
    if  centrality_kws is None:
        centrality_kws = {}

    centrality = c_func(p_graph, **centrality_kws)
    centrality = pd.Series(centrality).sort_index()
    if centrality.sum() == 0:
        LOG.error("type=unconnected_graph:frame={}. Try increasing the threshold.\n".format(frameidx+1))
        raise nx.exception.NetworkXException
    num_nodes = len(p_graph.nodes())
    nodes_axis = range(1, num_nodes + 1)
    if generate_plots:
        plt.plot(nodes_axis, centrality)
        plt.title("{0} {1}".format(prefix, label), fontsize=18)
        plt.xlabel("Node indices", fontsize=16)
        plt.ylabel("{}".format(label.capitalize()), fontsize=16)
        plt.tight_layout()
        plt.savefig("{}.png".format(prefix), dpi=300, bbox_inches='tight')
        plt.close()
    centrality = pd.Series(centrality)
    if label == "L":
        centrality.to_csv("{0}_avg_{1}.dat".format(prefix, label),
                          header=False, index=False)
    else:
        centrality.to_csv("{0}_{1}.dat".format(prefix, label), header=False, index=False)
    return centrality

def calc_drn_statistics(dataset, outfilename="DRN", drnstat="mean"):
    """Reporter/writer for DRN properties. Return a dictionary of DRN
    calculations mapped by their mean, median and std"""
    network_metrics = sorted(dataset.keys())
    drn_dict = {}
    drn_dict[drnstat] = pd.DataFrame(columns=network_metrics)
    for network_metric in network_metrics:
        data = dataset[network_metric]
        data = pd.DataFrame(data)
        if drnstat == "mean":
            drn_dict["mean"][network_metric] = data.mean(axis=1)
        elif drnstat == "median":
            drn_dict["median"][network_metric] = data.median(axis=1)
        elif drnstat == "std":
            drn_dict["std"][network_metric] = data.std(axis=1)
    drn_dict[drnstat].to_csv("{0}_{1}.csv".format(outfilename, drnstat),
                             index=False)
    return drn_dict

def topology_to_mmcif(infilename, outfilename, bfactors):
    """Converts a PDB file to mmcif with the DRN inserted as b-factors"""
    parser = PDBParser()
    pdb = parser.get_structure("topology", infilename)
    idx = 0
    for residue in pdb.get_residues():
        for atom in residue.get_atoms():
            atom.set_bfactor(bfactors[idx])
            idx += 1
    mmcifio = MMCIFIO()
    mmcifio.set_structure(pdb)
    mmcifio.save(outfilename)

def map_drn_to_structure(infilename, statistic=None, drn=None,
                         topology=None):
    """Maps mean, median or sd of DRN to a 3D structure """
    basename = os.path.splitext(infilename)[0]
    drn = drn[statistic]
    atom_filter = get_selection()
    topology = md.load_pdb(topology)
    topology.atom_slice(topology.top.select(atom_filter), inplace=True)
    topology.save_pdb("{}_top_.pdb".format(basename))
    for metric in drn.columns:
        outfilename = "{0}_{1}_{2}.cif".format(basename, statistic, metric)
        if drn[metric].sum() > 0:
            values = drn[metric].values
            topology_to_mmcif("{}_top_.pdb".format(basename), outfilename, values)
            LOG.info("Wrote DRN statistic as B-factor in {}\n".format(outfilename))
    os.remove("{}_top_.pdb".format(basename))

def calc_centralities(traj, traj_name, total_frames, args):
    """Calculates centralities across all MD frames and returns the summary
    statistics"""
    LOG.info("Calculating centralities...\n")
    dataset = {"BC":{}, "DC":{}, "CC":{}, "EC":{}, "ECC":{},
               "PR":{}, "katz":{}, "L":{}}
    atom_filter = get_selection(ligands=args.ligands)
    for frameidx, frame in enumerate(traj):
        frameidx += 1
        try:
            if total_frames:
                LOG.info("Progress: {}/{}\r".format(frameidx, total_frames))
            else:
                LOG.info("Progress: {} frames completed\r".format(frameidx))
            basename = ".".join(traj_name.split(".")[:-1])
            prefix = "{0}_{1}".format(basename, frame.time[0])
            p_graph = construct_graph(frame, atom_filter, prefix,
                                      args.threshold, args.save_graphs)
            if args.calc_BC:
                label = "BC"
                metric = calc_centrality(p_graph, c_func=nx.betweenness_centrality,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots,
                                         centrality_kws={"normalized":True})
                dataset[label][frameidx] = metric
            if args.calc_DC:
                label = "DC"
                metric = calc_centrality(p_graph, c_func=nx.degree_centrality,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
            if args.calc_CC:
                label = "CC"
                metric = calc_centrality(p_graph, c_func=nx.closeness_centrality,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
            if args.calc_EC:
                label = "EC"
                metric = calc_centrality(p_graph, c_func=nx.eigenvector_centrality_numpy,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
            if args.calc_ECC:
                label = "ECC"
                metric = calc_centrality(p_graph, c_func=nx.eccentricity,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
            if args.calc_PR:
                label = "PR"
                metric = calc_centrality(p_graph, c_func=nx.pagerank,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
            if args.calc_katz:
                label = "katz"
                metric = calc_centrality(p_graph, c_func=nx.katz_centrality,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
            if args.calc_L:
                label = "L"
                metric = calc_centrality(p_graph, c_func=calc_shortest_path,
                                         frameidx=frameidx, prefix=prefix, label=label,
                                         generate_plots=args.generate_plots)
                dataset[label][frameidx] = metric
        except nx.exception.NetworkXNoPath as nxnperror:
            sys.exit(nxnperror)
        except nx.exception.ExceededMaxIterations as nxemierror:
            LOG.error("type=max_iter_exceeded:frame={0}:exception={1}\n".format(frameidx, nxemierror))
            sys.exit(nxemierror)
        except Exception as ex:
            LOG.error("type=general:frame={0}:exception={1}\n".format(frameidx, Exception))
            sys.exit(ex)
    drn_dict = calc_drn_statistics(dataset, outfilename=basename, drnstat=args.drnstat)
    LOG.info("Saved DRN metrics as {}_*.csv.\n".format(basename))
    return drn_dict

def parse_args():
    """Argument parser"""
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
    parser.add_argument("--drnstat", default="mean", type=str,
                        help="The DRN statistic to compute [options: mean, median, std] \
                        (default=mean)")
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
                        help="Saves each frame's network graph in graphml and \
                              gml formats (default: disabled)")
    parser.add_argument("--lazy-load", action='store_true', default=False,
                        help="Read frames in a memory efficient way - use for \
                              big trajectories)")
    return parser

def main(args):
    """Program start"""
    if not any([args.calc_BC, args.calc_L, args.calc_DC, args.calc_ECC,
                args.calc_EC, args.calc_CC, args.calc_katz, args.calc_PR]):
        LOG.error("At least one of the network metrics (e.g. --calc-BC) must be set.\n")
        sys.exit(1)
    traj_name = os.path.basename(args.trajectory)
    traj, total_frames = load_trajectory(args.trajectory, args.topology,
                                         args.step, args.lazy_load)
    drn = calc_centralities(traj, traj_name, total_frames, args)
    map_drn_to_structure(os.path.splitext(traj_name)[0],
                         statistic=args.drnstat, drn=drn,
                         topology=args.topology)

LOG = Logger()

if __name__ == "__main__":
    CLI(parse_args(), main, LOG)
