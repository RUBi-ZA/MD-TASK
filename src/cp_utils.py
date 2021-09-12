#!/usr/bin/env python
#
# Additional methods for analysing coordination propensity data 
# 
#
# Script distributed under GNU GPL 3.0
#
# Author: David Penkler
# Date: 09-08-2018


import os
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib import cm
import matplotlib.pyplot as plt


class MDIterator:

    def __init__(self, traj_file, top, chunk=100, stride=1):
        self.iterator = md.iterload(traj_file, top=top,
                                    chunk=chunk, stride=stride)
        self.trajectory = None
        self.index = chunk - 1
        self.chunk = chunk

    def __iter__(self):
        return self

    def next(self):
        if self.index == self.chunk - 1:
            self.index = -1
            self.trajectory = self.iterator.next()

        self.index += 1

        try:
            return self.trajectory[self.index]
        except IndexError:
            raise StopIteration


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)
    

def plot_map(matrix, title, output_prefix, low=None, high=None,
             map_type='hot', axes=True, tick_interval=10):   
    '''Method for plotting an NxN heatmap of residue pair coordination propensities'''
    # Colour options
    #   inferno - black, purple, red, yellow  
    #   Blues = white to dark blue
    #   seismic - blue to red (reverse for red to blue)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["xtick.labelsize"] = 7
    plt.rcParams["ytick.labelsize"] = 7
    ax = sns.heatmap(matrix, cmap=map_type, xticklabels=tick_interval,
                     yticklabels=tick_interval)
    ax.set_xticklabels(ax.get_xticklabels(), fontdict={"rotation":90})
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Residue index")
    ax.invert_yaxis() 
    plt.savefig('%s.png' % output_prefix, bbox_inches="tight", dpi=300)
    plt.close('all')
