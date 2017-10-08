#!/usr/bin/env python
#
# Calculate the change in the valuses (BC or L) of each residue in a protein over the
# course of an MD simulation
#
# Script distributed under GNU GPL 3.0
#
# Author: David Brown
# Date: 17-11-2016

from natsort import natsorted

from lib.cli import CLI
from lib.utils import Logger
from lib.strategies import normalization

import numpy as np

import sys, argparse, matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def calc_delta(reference_file, alternative_files, normalizer, generate_plots=False):
    reference = np.loadtxt(reference_file)
    num_nodes = reference.shape[0]

    label = normalizer.get_label()

    alternatives = natsorted(alternative_files)

    log.info("Calculating %s for %d networks...\n" % (label, len(alternatives)))

    for i, alternative in enumerate(alternatives):
        log.info("Calculating %s (%d/%d)\r" % (label, i + 1, len(alternatives)))

        title = ".".join(alternative.split(".")[:-1])
        alternative = np.loadtxt(alternative)

        difference = alternative - reference
        difference = normalizer.normalize(difference, reference)
        prefix = "%s_%s_delta_%s" % (title, normalizer.get_prefix(), normalizer.matrix_type)

        np.savetxt("%s.dat" % prefix, difference)

        if generate_plots:
            node_axis = range(1, num_nodes + 1)
            plt.plot(node_axis, difference)
            plt.axhline(0, color='black')
            plt.title("%s %s" % (title, label), fontsize=18)
            plt.xlabel('Residue Numbers', fontsize=16)
            plt.ylabel(label, fontsize=16)
            plt.savefig("%s.png" % prefix, dpi=300, bbox_inches="tight")
            plt.close()

    log.info("\n")


def get_normalizer(matrix_type, mode):
    type_map = {
        "L": "standard",
        "BC": "plusone"
    }

    if matrix_type not in type_map:
        log.info("ERROR: --matrix-type must be specified and must be either 'BC' or 'L'")
        sys.exit(1)

    normalization_mode = mode if mode else type_map[matrix_type]

    return getattr(normalization, normalization_mode, "standard")(matrix_type)


def main(args):
    if args.normalize:
        normalizer = get_normalizer(args.matrix_type, args.normalization_mode)
    else:
        normalizer = normalization.none(args.matrix_type)

    calc_delta(args.reference, args.alternatives, normalizer, args.generate_plots)


log = Logger()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--matrix-type", help="The type of values in the matrices i.e. BC or L", default=None)
    parser.add_argument("--reference", help="The reference matrix (.dat)")
    parser.add_argument("--alternatives", help="The alternative matrices (.dat)", nargs="*")
    parser.add_argument("--normalize", help="Normalizes the values", action='store_true', default=False)
    parser.add_argument('--normalization-mode', help="Method used to normalize (default for L = standard, default for BC = plusone)", default=None)
    parser.add_argument("--generate-plots", help="Plot results - without setting this flag, no graph will be generated", action='store_true', default=False)

    CLI(parser, main, log)
