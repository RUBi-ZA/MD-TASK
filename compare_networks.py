#!/usr/bin/env python
#
# Compare network measurements such as BC and L by plotting them on a graph
#
# Script distributed under GNU GPL 3.0
# 
# Author: David Brown
# Date: 17-11-2016

import numpy as np
import matplotlib.pyplot as plt

import sys, argparse

    

def plot_comparison(reference, alternative, reference_label, alternative_label, prefix="compare", y_label=None, ylim=None):
    difference = alternative - reference
    title = "%s vs %s ($\Delta$%s)" % (reference_label, alternative_label, y_label)
    
    num_nodes = reference.shape[0]
    node_axis = range(1, num_nodes + 1)
    
    plt.subplots(figsize=(15, 20))
    
    #comparison plot
    plt.subplot(211)
    plt.axhline(0, color='black')
    
    lines = []
    
    lines.append(plt.plot(node_axis, reference, "black", label=reference_label)[0])
    lines.append(plt.plot(node_axis, alternative, "red", label=alternative_label)[0])
    
    plt.legend(handles=lines)
    plt.title("%s vs %s (%s)" % (alternative_label, reference_label, y_label), fontsize=18)
    plt.xlabel('Residue Numbers', fontsize=16)
    plt.ylabel(y_label, fontsize=16)

    #delta plot
    plt.subplot(212)
    plt.axhline(0, color='black')
    
    line, = plt.plot(node_axis, difference, "red", label="$\Delta$%s" % y_label)
    plt.legend(handles=[line])
    plt.title("%s - %s ($\Delta$%s)" % (alternative_label, reference_label, y_label), fontsize=18)
    plt.xlabel('Residue Numbers', fontsize=16)
    plt.ylabel("$\Delta$%s" % y_label, fontsize=16)
    
    if ylim:
        axes = plt.gca()
        axes.set_ylim(ylim)
    
    filename = "%s_comp.png" % prefix
    log("Generate plot: %s" % filename)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


def main(args):
    reference = np.loadtxt(args.reference)
    alternative = np.loadtxt(args.alternative)
    
    y_label = args.y_label
    if not y_label:
        y_label = "Measurement"
    
    ylim = None
    if args.y_max is not None and args.y_min is not None:
        ylim = [float(args.y_min), float(args.y_max)]
    
    plot_comparison(reference, alternative, args.reference_label, args.alternative_label, args.prefix, y_label, ylim)


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
    parser.add_argument("--reference", help="The reference network (.dat)")
    parser.add_argument("--alternative", help="The alternative network (.dat)")
    parser.add_argument("--prefix", help="Prefix for output files")
    parser.add_argument("--reference-label", help="Label to display on graphs for reference network", default="")
    parser.add_argument("--alternative-label", help="Label to display on graphs for alternative network", default="")
    parser.add_argument("--y-label", help="Label for y-axis", default=None)
    parser.add_argument("--y-max", help="Maximum value on y-axis", default=None)
    parser.add_argument("--y-min", help="Minimum value on y-axis", default=None)
    
    args = parser.parse_args()
    
    #set up logging
    silent = args.silent
    
    if args.log_file:
        stream = open(args.log_file, 'w')
    
    #run script
    main(args)
    
    #close logging stream
    stream.close()
    
 
