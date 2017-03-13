#!/bin/env python

import sys
import gzip
import argparse
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

def set_from_string(string):
    '''
    Convert a comma separated string to a set of strings
    '''
    if string == '':
        return set()
    else:
        return set(string.split(","))

def overlapping_ids(query_set, filtering_sets):
    '''
    Return the ids for all sets in our list of possible overlapping sets.
    filtering_sets is expected to be a tuple where index 0 is the name and index 1 is the set
    '''
    found = []
    for f in filtering_sets:
        if query_set & f[1]:
            found.append(f[0])
    return found

def load_filter_file(filter_file):
    '''
    Read the file we're going to use as a filter to determine if lines should be output.
    This returns a list containing tuples where the first item is the variant id and the second is the set of ids from sname.
    '''
    filter_list = list()

    vcf = Vcf()
    header_lines = list()
    in_header = True
    for line in filter_file:
        if in_header:
            header_lines.append(line)
            if line[0:6] == '#CHROM':
                in_header = False
                vcf.add_header(header_lines)
        else:
            v = line.rstrip().split('\t')
            var = Variant(v, vcf)
            filter_list.append((var.var_id, set_from_string(var.get_info('SNAME'))))
    return filter_list

def sname_filter(input_stream, filter_file, output_stream, complement):
    '''
    This reads a VCF stream, determines if the line overlaps any from the filter_file by sname and outputs.
    '''
    filter_list = load_filter_file(filter_file)

    vcf = Vcf()
    in_header = True    
    header_lines = list()
    sample_list = None
    for line in input_stream:
        if in_header:
            header_lines.append(line)
            if line[0:6] == '#CHROM':
                in_header = False
                vcf.add_header(header_lines)
                vcf.add_info('FOUND', '.', 'String', 'Variant id in other file')
                output_stream.write(vcf.get_header() + '\n')
        else:
            v = Variant(line.rstrip().split('\t'), vcf)
            sname_set = set_from_string(v.get_info('SNAME'))
            found = overlapping_ids(sname_set, filter_list)
            if bool(found) != complement:
                v.set_info('FOUND', ','.join(found))
                output_stream.write(v.get_var_string() + '\n')

def description():
    return 'look for variants sharing the same original call between two VCF files'

def add_arguments_to_parser(parser):
    parser.add_argument('-v', '--complement', action='store_true', dest='complement', default=False, help='return complement of overlap')
    parser.add_argument("-i", "--input", dest="input", metavar='<VCF>', help="VCF file containing variants to be output")
    parser.add_argument("-f", "--filter", dest="filter_file", metavar='<VCF>', type=argparse.FileType('r'), help="VCF file containing variants used to determine if site should be output")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<VCF>', default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        return sname_filter(stream, args.filter_file, args.output, args.complement)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
