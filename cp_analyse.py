# Methods for analysing coordination propensity data 
# 
#
# Script distributed under GNU GPL 3.0
#
# Author: David Penkler
# Date: 09-08-2018

import os, sys, argparse

import numpy as np

import cp_utils
from cp_utils import *
from lib.cli import CLI
from lib.utils import Logger

def calculate_avg_cp(matrix, prefix):
    '''Calcualtes the average cp for each residue, 
       where average cp is the average cp between residue i 
       and all other residues'''
       
    dimensions = len(matrix[:,0])
    ave = []
    for res in range(0, dimensions):
    	ave.append(np.mean(matrix[res,:]))
    	
    np.savetxt('%s_avg_cp.csv' % prefix, ave, delimiter=",")

def calculate_sliding_avg(matrix,prefix,window=4):
    '''Calculates the average cp for each residue using a sliding window,
       where average cp is the average cp between residue i and neigbouring 
       residues +- window away'''
       
    dimensions = len(matrix[:,0])

    total = 0
    count = 0
    totals = []
    
    for row in range(0,dimensions):
        if row >= window and row <= dimensions-window:
            avg = np.mean(matrix[row,row-window:row+window])
            total += avg
            totals.append(avg)
            
        elif row < window:
            avg = np.mean(matrix[row,0:row+window])
            total += avg
            totals.append(avg)

        elif row > dimensions-window:
            avg = np.mean(matrix[row,row-window:dimensions])
            total += avg
            totals.append(avg)
   
    np.savetxt('%s_sliding_window_avg_cp.csv' % prefix, totals, delimiter=',')
    
                
def calculate_difference(m1,m2,prefix):
    '''Calculates the difference matrix and plots it as a heatmap
       diff = m1 - m2'''
    diff = m1-m2
    np.save('%s.npy'%(prefix),diff)
    
    low = np.min(diff)
    high = np.max(diff)
    
    cp_utils.plot_map(diff,'Delta CP','%s_delta_cp_heatmap' % args.prefix,low=low, high=high, map_type='seismic', axes=False)
    return diff
  


def main(args):
    log.info(" -Loading coordination propensity data ...")
    
    cp_matrix_1 = np.load(args.cp)
    
    if args.avg:
        log.info(" -Calculating average CP for each residue ...")
        calculate_avg_cp(cp_matrix_1, args.prefix) 
    
    if args.sliding_avg:
        log.info(" -Calculating average CP for each residue using a sliding window of %s ..." % args.window)
        calculate_sliding_avg(cp_matrix_1,args.prefix,window=args.window)
    
    if args.diff:
        log.info(" -Calculating and plotting delta CP heatmap...")
        cp_matrix_2 = np.load(args.diff)
        diff_matrix = calculate_difference(cp_matrix_1, cp_matrix_2, args.prefix)

log = Logger()

if __name__ == "__main__":

    #parse cmd arguments
    parser = argparse.ArgumentParser()

    #custom arguments
    parser.add_argument("cp", help="Coordination propensity matrix (matrix.npy)")
    parser.add_argument("--avg", help="Calculate the total average cp for each residue", action="store_true")
    parser.add_argument("--sliding_avg", help="Calculate average cp for every residue using a sliding window", action='store_true')
    parser.add_argument("--window", help="Sliding window", type=int, default=4)
    parser.add_argument("--diff", help="Calculate the difference between cp and diff matricies", default=None)
    
    parser.add_argument("--prefix", help="Prefix for output files (default: cp)", default="cp")

    #args = parser.parse_args()
	
    #run script
    #main(args)

    CLI(parser, main, log)
	