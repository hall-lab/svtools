import sys
import argparse
from svtools.external_cmd import ExternalCmd

class VcfSort(ExternalCmd):
    def __init__(self):
        super(VcfSort, self).__init__('vcfsort', 'bin/vcfsort')

def description():
    return 'sort a VCF file'

def add_arguments_to_parser(parser):
    parser.add_argument('input', metavar='<VCF>', nargs='?', help='VCF file to sort (default: stdin)')
    parser.add_argument('output', metavar='<VCF>', nargs='?', help='output file to write to (default: stdout)')
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

    sort_cmd_runner = VcfSort()
    sort_cmd_runner.run_cmd_with_options(opts)

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
