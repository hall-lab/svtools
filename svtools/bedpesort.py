import sys
import argparse
from svtools.external_cmd import ExternalCmd

class BedpeSort(ExternalCmd):
    def __init__(self):
        super(BedpeSort, self).__init__('bedpesort', 'bin/bedpesort')

def description():
    return 'sort a BEDPE file'

def add_arguments_to_parser(parser):
    parser.add_argument('input', metavar='<BEDPE file>', nargs='?', help='BEDPE file to sort')
    parser.add_argument('output', metavar='<output file>', nargs='?', help='output file to write to')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    opts = list()
    if args.input:
        opts.append(args.input)
    if args.output:
        opts.append(args.output)

    sort_cmd_runner = BedpeSort()
    sort_cmd_runner.run_cmd_with_options(opts)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
