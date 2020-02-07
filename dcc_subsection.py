#Caroline Ross 14 December 2018
#Plots a sub-section of dcc correlation
#Input = the correlation.txt file from MD-TASK

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import numpy as np


def plot_map(correlation, title, output_prefix, x_labels, y_labels):
    M = np.array(correlation)

    fig, ax = plt.subplots()
    colors = [('white')] + [(cm.jet(i)) for i in range(40,250)]

    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(M, cmap=new_map, vmin=-1, vmax=1)

    fig = plt.gcf()
    ax.set_frame_on(False)
    ax.grid(False)
    
    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(M.shape[0])+0.5, minor=False)
    ax.set_xticks(np.arange(M.shape[1])+0.5, minor=False)

    ax.set_xticklabels(x_labels,fontsize=8, minor=False) 
    ax.set_yticklabels(y_labels,fontsize=8, minor=False)

    plt.xticks(rotation=90)

    # Turn off all the ticks
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
 

    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    
    plt.title(title, fontsize=16)
    plt.xlabel('Residue Index', fontsize=12)
    plt.ylabel("Residue Index", fontsize=12)

    cbar = plt.colorbar(heatmap, orientation="vertical")
    plt.savefig('%s.png' % output_prefix, dpi=300)
    plt.close('all')


def print_correlation(correlation, output_prefix):
    with open("%s.txt" % output_prefix, "w") as w:
        rows = correlation.shape[0]
        cols = correlation.shape[1]

        for r in range(rows):
            for c in range(cols):
                w.write('%s ' % str(correlation[r,c]))
            w.write('\n')

#Reads in correlation.txt file:
try:
    f = open('correlation.txt', 'r') #change file name here for use on different file
    correlation_values = f.readlines() #reads in all lines of correlation.txt
    f.close() #close file
except IOError:
    print ('\n**************************************\nERROR!! FILE NOT FOUND:\n**************************************\n') #error if correlation.txt file not found
    sys.exit()


#############################################################################################################################################################
# This code can be updated to analyse any section of the protein
#
#############################################################################################################################################################
#Extract a section

proteinSectionA_start = 1 #change as required
proteinSectionA_end = 27

proteinSectionB_start = 109
proteinSectionB_end = 120

xatoms = proteinSectionA_end-proteinSectionA_start+1 #(sectionA will be rows in matrix - labeled along y axis)
yatoms = proteinSectionB_end-proteinSectionB_start+1 #(sectionB will be columns in matrix - labeled along x axis)

Sub_cMatrix = np.zeros((xatoms, yatoms)) #set size of sub_matrix

for i,x in enumerate(range(proteinSectionA_start-1,proteinSectionA_end)):
    atom_specific_correlation = correlation_values[x].split()
    for j, y in enumerate(range(proteinSectionB_start-1,proteinSectionB_end)):
       x_yCorrelation = float(atom_specific_correlation[y].strip())
       Sub_cMatrix[i, j] = x_yCorrelation

x_labels = []
for i in range(proteinSectionB_start,proteinSectionB_end+1):
    x_labels.append(str(i))

y_labels = []
for i in range(proteinSectionA_start,proteinSectionA_end+1):
    y_labels.append(str(i))

Title = "Sub-Correlation Plot" #Change here as required
ImageName = "Sub-correlation" #Change here as required

plot_map(Sub_cMatrix, Title,ImageName, x_labels, y_labels)
