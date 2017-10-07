#!/usr/bin/env python
#
# Calculate the change in the average shotest paths of each residue over the
# course of an MD simulation
#
# Script distributed under GNU GPL 3.0
#
# Author: David Brown
# Date: 17-11-2016

from datetime import datetime
import sys, argparse, calc_delta

from lib.utils import format_seconds


def main(args):
    args.matrix_type = "L"
    calc_delta.main(args)


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
    parser.add_argument('--normalization-mode', help="Method used to normalize (default: (Delta L/L))", default=None)
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

