#!/usr/bin/env python
#
# Calculate average network measurements over the course of the simulation as
# well as standard deviation of the measurements
#
# Script distributed under GNU GPL 3.0
# 
# Author: David Brown
# Date: 17-11-2016

import numpy as np

from datetime import datetime

from lib.utils import *

import os, sys, argparse, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def split_array(arr, pos):
    return arr[:pos], arr[pos:]
    

def combine_arrays(array_files):
    matrix = None
    for network in array_files:
        if matrix is None:
            matrix = np.loadtxt(network)
        else:
            matrix = np.vstack([matrix, np.loadtxt(network)])
    
    return matrix


def plot_graph(network, err=None, start_x=1, color="black", ecolor="red", title="Title", x_label="X", y_label="Y", ylim=None):
    start_x = int(start_x)
    
    num_nodes = network.shape[0]
    nodes_axis = range(start_x, num_nodes + start_x)
    
    plt.axhline(0, color='black')
    
    if err is not None:
        plt.errorbar(nodes_axis, network, err, color="black", ecolor="red")
    else:
        plt.plot(nodes_axis, network, color="black")
    
    if ylim:
        axes = plt.gca()
        axes.set_ylim(ylim)
    
    plt.title(title, fontsize=18)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)


def plot_graph_split(split_pos, avg_matrix, std, initial_x_1, initial_x_2, title_1, title_2, x_label, y_label, ylim=None):
    num_nodes = avg_matrix.shape[0]
    pos = int(split_pos)
    
    avg_matrix_1,avg_matrix_2 = split_array(avg_matrix, pos)
    std_1,std_2 = split_array(std, pos)
    
    plt.subplots(figsize=(25, 20))
    
    #first plot
    plt.subplot(211)
    
    plot_graph(avg_matrix_1, std_1, args.initial_x_1, title=title_1,
        x_label=x_label, y_label=y_label, ylim=ylim)
    
    #second plot
    plt.subplot(212)
    
    plot_graph(avg_matrix_2, std_2, args.initial_x_2, title=title_2,
        x_label=x_label, y_label=y_label, ylim=ylim)


def main(args):
    prefix = args.prefix
    start_index = args.initial_x
    
    if not args.data_type in ["BC", "delta-BC", "L", "delta-L"]:
        log("Unrecognized data type. Exiting...\n")
        sys.exit(1)
    
    # average BC
    if args.data_type == "BC":
        prefix += "_BC"
    
    # delta average BC    
    elif args.data_type == "delta-BC":
        prefix += "_delta_BC"
        
    # average L
    elif args.data_type == "L":
        prefix += "delta_L"
        
    # delta average L
    elif args.data_type == "delta-L":
        prefix += "_delta_L"
    
    # calculate matrices    
    matrix = combine_arrays(args.data)
    np.savetxt("%s_combined.dat" % prefix, matrix)
    
    std = matrix.std(0)
    np.savetxt("%s_std_dev.dat" % prefix, std)
    
    avg_matrix = np.mean(matrix, axis=0)
    np.savetxt("%s_avg.dat" % prefix, avg_matrix)
    
    if args.generate_plots:
        log("Generating graph: %s_avg.png\n" % prefix)
        
        ylim = None
        if args.y_max is not None and args.y_min is not None:
            ylim = [float(args.y_min), float(args.y_max)]
        
        if args.split_pos:
            plot_graph_split(
                int(args.split_pos), avg_matrix, std, args.initial_x_1, args.initial_x_2, 
                args.title_1, args.title_2, args.x_label, args.y_label, ylim=ylim
            )
            
        else:  
            plot_graph(
                avg_matrix, std, args.initial_x, title=args.title, x_label=args.x_label,
                y_label=args.y_label, ylim=ylim
            )
        
        plt.savefig("%s_avg.png" % prefix, dpi=300, bbox_inches='tight')
        plt.close()


silent = False
stream = sys.stdout

def log(message):
    global silent
    global stream
    
    if not silent:
        stream.write(message)


if __name__ == "__main__":
    
    #parse cmd arguments
    parser = argparse.ArgumentParser()
    
    #standard arguments for logging
    parser.add_argument("--silent", help="Normalizes the values (Delta L/L)", action='store_true', default=False)
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    
    #custom arguments
    parser.add_argument("--data", help="The .dat files that will be averaged", nargs="*")
    parser.add_argument("--data-type", help="Type of data - BC/delta-BC/L/delta-L")
    
    parser.add_argument("--prefix", help="Prefix used to name outputs", default="network")
    
    parser.add_argument("--generate-plots", help="Generate figures/plots", action='store_true', default=False)
    
    # plot arguments (only used with --generate-plots)
    parser.add_argument("--x-label", help="Label for x-axis (use $\Delta$ for delta sign)", default=None)
    parser.add_argument("--y-label", help="Label for y-axis (use $\Delta$ for delta sign)", default=None)
    parser.add_argument("--y-max", help="Maximum value on y-axis", default=None)
    parser.add_argument("--y-min", help="Minimum value on y-axis", default=None)
    
    # argument if generating a single plot
    parser.add_argument("--title", help="Title of plot (use $\Delta$ for delta sign)", default=None)
    parser.add_argument("--initial-x", help="The start index of the x-axis", default=1)
    
    # arguments if splitting the plot in two
    parser.add_argument("--split-pos", help="Position to split the network at (for large networks)", default=None)
    parser.add_argument("--title-1", help="Title of first network", default="Plot 1")
    parser.add_argument("--title-2", help="Title of second network", default="Plot 2")
    parser.add_argument("--initial-x-1", help="The start index of the x-axis for the first network", default=1)
    parser.add_argument("--initial-x-2", help="The start index of the x-axis for the second network", default=1)
    
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
    
    log("Completed at: %s\n" % str(end))
    log("- Total time: %s\n" % str(time_taken))
    
    #close logging stream
    stream.close()
