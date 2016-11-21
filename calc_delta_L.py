#!/usr/bin/env python
#
# Calculate the change in the average shotest paths of each residue over the 
# course of an MD simulation
#
# Script distributed under GNU GPL 3.0
# 
# Author: David Brown
# Date: 17-11-2016


from natsort import natsorted
from datetime import datetime

from lib.utils import *

import numpy as np
import matplotlib.pyplot as plt

import os, sys, argparse



def calc_delta_L(reference_file, alternative_files, normalize=False, generate_plots=False):
    reference = np.loadtxt(reference_file)
    num_nodes = reference.shape[0]
    
    if normalize:
        label = "L/L"
    else:
        label = "L"
    
    delta_L = np.zeros((num_nodes, num_nodes))
    alternatives = natsorted(alternative_files)
    
    log("Calculating delta %s for %d networks..." % (label, len(alternatives)))
    
    for i, alternative in enumerate(alternatives):
        log("Calculating delta %s between %s and %s (%d/%d)" % (label, alternative, reference_file, i + 1, len(alternatives)))
        
        prefix = ".".join(alternative.split(".")[:-1])
        alternative = np.loadtxt(alternative)
        
        difference = alternative - reference
        if normalize:
            difference = difference / reference
            prefix += "_norm"
        
        delta_L[i,:] = difference
        
        np.savetxt("%s_delta_L.dat" % prefix, difference)
        
        if generate_plots:
            node_axis = range(1, num_nodes + 1)
            plt.plot(node_axis, difference)
            plt.axhline(0, color='black')
            plt.title("%s $\Delta$ %s" % (prefix, label), fontsize=18)
            plt.xlabel('Residue Numbers', fontsize=16)
            plt.ylabel("$\Delta$ %s" % label, fontsize=16)
            plt.savefig("%s_delta_L.png" % prefix, dpi=300, bbox_inches='tight')
            plt.close()
            
            log("Plot generated: %s_delta_L_avg.png" % prefix)
    
    log("Completed.")
    
    return delta_L



def main(args):
    calc_delta_L(args.reference, args.alternatives, args.normalize, args.generate_plots)



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
    parser.add_argument("--reference", help="The reference avg L matrix (.dat)")
    parser.add_argument("--alternatives", help="The alternative avg L matrices (.dat)", nargs="*")
    parser.add_argument("--normalize", help="Normalizes the values (Delta L/L)", action='store_true', default=False)
    parser.add_argument("--generate-plots", help="Plot results - without setting this flag, no graph will be generated", action='store_true', default=False)
    
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
    