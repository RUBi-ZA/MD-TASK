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
from lib.strategies import normalization

import numpy as np

import os, sys, argparse, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt



def calc_delta_L(reference_file, alternative_files, normalize, normalization_function, generate_plots=False):
    reference = np.loadtxt(reference_file)
    num_nodes = reference.shape[0]
    
    if normalize:
        label = "L/L"
    else:
        label = "L"
    
    alternatives = natsorted(alternative_files)
    
    log("Calculating delta %s for %d networks...\n" % (label, len(alternatives)))
    
    for i, alternative in enumerate(alternatives):
        log("Calculating delta %s (%d/%d)\r" % (label, i + 1, len(alternatives)))
        
        prefix = ".".join(alternative.split(".")[:-1])
        alternative = np.loadtxt(alternative)
        
        difference = alternative - reference
        if normalize:
            difference = normalization_function(difference, reference)
            prefix += "_norm"
        
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
    
    log("\n")



def main(args):
    normalization_function = getattr(normalization, args.normalization_mode, 'default')
    calc_delta_L(args.reference, args.alternatives, args.normalize, normalization_function, args.generate_plots)



silent = False
stream = sys.stdout

def log(message):
    global silent
    global stream
    
    if not silent:
        stream.write(message)
        stream.flush()



if __name__ == "__main__":
    
    #parse cmd arguments
    parser = argparse.ArgumentParser()
    
    #standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    
    #custom arguments
    parser.add_argument("--reference", help="The reference avg L matrix (.dat)")
    parser.add_argument("--alternatives", help="The alternative avg L matrices (.dat)", nargs="*")
    parser.add_argument("--normalize", help="Normalizes the values", action='store_true', default=False)
    parser.add_argument('--normalization-mode', help="Method used to normalize (default: (Delta L/L))", nargs='?', const="default", default="default")
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
    
    log("Completed at: %s\n" % str(end))
    log("- Total time: %s\n" % str(time_taken))
    
    #close logging stream
    stream.close()
    
