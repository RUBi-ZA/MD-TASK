#!/usr/bin/env python3
#
# Calculate the change in the average shotest paths of each residue over the
# course of an MD simulation
#
# Script distributed under GNU GPL 3.0
#
# Author: David Brown
# Date: 17-11-2016

import argparse, calc_delta

from lib.cli import CLI
from lib.utils import Logger


def main(args):
    args.matrix_type = "L"
    calc_delta.main(args)


log = Logger()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--reference", help="The reference avg L matrix (.dat)")
    parser.add_argument("--alternatives", help="The alternative avg L matrices (.dat)", nargs="*")
    parser.add_argument("--normalize", help="Normalizes the values", action='store_true', default=False)
    parser.add_argument('--normalization-mode', help="Method used to normalize (default: (Delta L/L))", default=None)
    parser.add_argument("--generate-plots", help="Plot results - without setting this flag, no graph will be generated", action='store_true', default=False)

    CLI(parser, main, log)
