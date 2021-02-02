#!/usr/bin/env python
"""
Compare multiple weighted contact heatmap around a common residue position
"""
__author__ = "Olivier Sheik Amamuddy"
__copyright__ = "Copyright 2020, Research Unit in Bioinformatics"
__license__ = "GNU GPL 3.0"
__version__ = "1.3.2"
__email__ = "oliserand@gmail.com"
__status__ = "Production"
__date__ = "15th December 2020"

import os
import sys
import argparse
from datetime import datetime
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
from lib.utils import format_seconds

def plot_heatmap(allnetworks_copy, outfilename, xlabel='Target residue',
                ylabel='Contacting residues', figwidth=16,
                figheight=16, fontsize=12, dpi=450, annotate=True,
                shorten_yticks=False, xtickfontsize=6, ytickfontsize=6):
    """Plots the heat map"""
    plt.figure(figsize=[figwidth, figheight])
    axes = sns.heatmap(data=allnetworks_copy.transpose(), cmap='hot_r',
                       annot=annotate, linecolor="lightgrey",
                       linewidth=0.5, xticklabels=1, yticklabels=1,
                       annot_kws={"rotation":90})
    axes.set_yticklabels(axes.get_yticklabels(), rotation=0, fontsize=ytickfontsize)
    axes.set_xticklabels(axes.get_xticklabels(), rotation=90, fontsize=xtickfontsize)
    axes.set_ylabel(xlabel, fontsize=fontsize)
    axes.set_xlabel(ylabel, fontsize=fontsize)
    if shorten_yticks:
        new_labels = [i.get_text().split()[-1] for i in axes.get_yticklabels()]
        axes.set_yticklabels(new_labels)
    plt.tight_layout()
    plt.savefig(outfilename, dpi=dpi, bbox_inches="tight")

def main(arg):
    """Script main"""
    networks        = arg.networks
    outfilename     = arg.output
    allnetworks     = {}
    common_residues = set()
    plt.rcParams["font.family"] = "serif"
    for networkfile in networks:
        try:
            variant = pd.read_csv(networkfile, header=None,
                    names=['source','dest','weight'])
        except FileNotFoundError as _:
            log("\nERROR: {} file not found.\n".format(networkfile))
            sys.exit(1)
        center = variant.source.values[0]
        label_prefix = os.path.basename(networkfile).replace("_network.csv","").replace(".csv","")
        center = "({}) {}".format(label_prefix, center)
        allnetworks[center] = variant
        _ = [common_residues.add(residue) for residue in variant.dest.values]
    common_residue_order = pd.Series(
            [int(i[3:-2]) for i in common_residues]).sort_values(ascending=True).index
    common_residues = pd.Series(list(common_residues)).iloc[common_residue_order].values
    allnetworks_copy = pd.DataFrame(data=0,columns=allnetworks.keys(),
                                    index=common_residues)
    for center in allnetworks:
        for common_residue in common_residues:
            if common_residue in allnetworks[center].dest.values:
                allnetworks_copy.loc[common_residue, center] = \
                        allnetworks[center].loc[
                                allnetworks[center].dest == common_residue, 'weight'].values
    plot_heatmap(allnetworks_copy, outfilename, figheight=arg.figheight,
                 figwidth=arg.figwidth, fontsize=arg.fontsize, dpi=arg.dpi,
                 annotate=args.annotate, shorten_yticks=args.shorten_yticks,
                 xtickfontsize=args.xtickfontsize, ytickfontsize=args.ytickfontsize)

def check_filename(filename):
    """Sanity check for filenames"""
    if " " in filename:
        return False
    return True

def parse_args():
    """parse arguments"""
    parser = argparse.ArgumentParser("Heat map representation for comparing"\
            "several contact maps around a common residue position in "\
            "multiple mutants.")
    #standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging",
                        action='store_true', default=False)
    parser.add_argument("--log-file", default=None,
                        help="Output log file (default: standard output)")
    #custom arguments
    parser.add_argument('networks', nargs='*',
                        help="The network CSV files generated from running\
                        the 'contact_map.py' tool")
    parser.add_argument("--output", default="contact_heatmap.png",
                        help="Output PNG filename (default:contact_heatmap.png)")
    parser.add_argument("--annotate", action="store_true",
                        help="Add numerical annotations to heat map (default:no)")
    parser.add_argument("--dpi", default=450, type=int,
                        help="Resolution (in DPI) for the PNG file\
                        (default:450)")
    parser.add_argument("--figwidth", default=8, type=int,
                        help="Figure width (in inches)\
                        (default:8)")
    parser.add_argument("--figheight", default=5, type=int,
                        help="Figure height (in inches)\
                        (default:5)")
    parser.add_argument("--fontsize", default=12, type=int,
                        help="Font size (in points)\
                        (default:12)")
    parser.add_argument("--xtickfontsize", default=4.5, type=int,
                        help="Font size for x ticks (in points)\
                        (default:5)")
    parser.add_argument("--ytickfontsize", default=5, type=int,
                        help="Font size for y ticks (in points)\
                        (default:5)")
    parser.add_argument("--shorten_yticks", action="store_true",
                        help="Shorten y tick labels by ignoring file names\
                        (default:no)")
    arguments = parser.parse_args()
    return arguments

SILENT = False
STREAM = sys.stdout

def log(message):
    """The logger"""
    global SILENT
    global STREAM
    if not SILENT:
        STREAM.write(message)

if __name__ == "__main__":
    args = parse_args()
    #set up logging
    SILENT = args.silent
    if args.log_file:
        STREAM = open(args.log_file, 'w')
    if len(args.networks) < 2:
        log("ERROR: At least two CSV contact files generated from contact_map.py are "\
            "required. Try the -h option.\n")
        sys.exit(1)
    if not all(map(check_filename, args.networks)):
        log("ERROR: Make sure there are no spaces in the network CSV file names.\n")
        sys.exit(1)
    start = datetime.now()
    log("INFO: Started at: %s\n" % str(start))
    #run script
    main(args)
    end = datetime.now()
    time_taken = format_seconds((end - start).seconds)
    log("INFO: Completed at: %s\n" % str(end))
    log("INFO: - Total time: %s\n" % str(time_taken))
    #close logging stream
    STREAM.close()
