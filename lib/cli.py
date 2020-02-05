from datetime import datetime
from utils import format_seconds

def CLI(parser, callback, log):
    #standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)

    args = parser.parse_args()

    #set up logging
    log.silent = args.silent
    if args.log_file:
        log.stream = open(args.log_file, 'w')

    start = datetime.now()
    log.info("Started at: %s\n" % str(start))

    #run script
    callback(args)

    end = datetime.now()
    time_taken = format_seconds((end - start).seconds)

    log.info("Completed at: %s\n" % str(end))
    log.info("- Total time: %s\n" % str(time_taken))
