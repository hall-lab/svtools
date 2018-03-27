#!/usr/bin/env python

from __future__ import division
import argparse
import sys
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

class VCFReader(object):
    def __init__(self, stream):
        self.vcf_obj = Vcf()
        self.stream = stream
        header = list()
        for line in stream:
            if line[0] != '#':
                raise RuntimeError('Error parsing VCF header. Line is not a header line. {}'.format(line))
            header.append(line)
            if line.startswith('#CHROM\t'):
                # end of header
                break
        self.vcf_obj.add_header(header)

    def __iter__(self):
        for line in self.stream:
            yield Variant(line.rstrip().split('\t'), self.vcf_obj)


def load_deletion_sizes(stream):
    minimum_del_size = dict()
    for line in stream:
        if line.startswith('Sample\t'):
            continue
        sample, size, overlap = line.rstrip().split('\t')
        if sample not in minimum_del_size:
            minimum_del_size[sample] = abs(int(size))
        else:
            raise RuntimeError('Size for {0} already set. Does your file of sizes include multiple lines with the same sample name?'.format(sample))
    return minimum_del_size, max(minimum_del_size.values())

def set_missing(input_stream, deletion_sizes, output_stream, max_del_size, sr_cutoff):
    valid_types = set(('DEL', 'MEI'))
    for variant in input_stream:
        if variant.get_info('SVTYPE') in valid_types:
            # NOTE this will raise an exception if SVLEN is null
            length = abs(int(variant.get_info('SVLEN')))
            if length < max_del_size:
                split_read_support = 0
                total_depth = 0
                for s in variant.sample_list:
                    g = variant.genotype(s)
                    if g.get_format('GT') not in ('./.', '0/0'):
                        split_read_support += int(g.get_format('AS'))
                        total_depth += int(g.get_format('DP'))
                if total_depth > 0 and (split_read_support / total_depth) < sr_cutoff:
                    # Only set to null if PE support is our only source of
                    # information. Not counting soft-clips here.
                    # This may be a bad idea as even with SR support the
                    # lack of power to detect PE reads could skew us away from 
                    # the correct genotype.
                    # A better method might be to regenotype using only 
                    # split-read support if the SV is too small.
                    logged = False
                    for sample in variant.sample_list:
                        if sample in deletion_sizes and length < deletion_sizes[sample]:
                            gt = variant.genotype(sample)
                            gt.set_format('GT', './.')
                            if not logged:
                                sys.stderr.write('Applying small deletion filter to {0}\n'.format(variant.var_id))
                                logged=True

        output_stream.write(variant.get_var_string())
        output_stream.write('\n')
                
def description():
    return 'set genotypes of deletions smaller than a per-sample cutoff to missing if no splitread support in the sample'

def add_arguments_to_parser(parser):
    parser.add_argument("-i", "--input", required=True, dest="input", metavar='<VCF>', help="VCF file containing variants to be output")
    parser.add_argument("-t", "--thresholds", required=True, dest="threshold_file", metavar='<TXT>', type=argparse.FileType('r'), help="Tab-separated file of sample name and minimum deletion size used to determine if site should be output")
    parser.add_argument("-s", "--split-read-fraction", required=True, dest="sr_cutoff", metavar='<FLOAT>', type=float, help="Minimum fraction of split read support for the site to be excluded from filtering.")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<VCF>', default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    deletion_size_map, max_size = load_deletion_sizes(args.threshold_file)
    with su.InputStream(args.input) as stream:
        variant_stream = VCFReader(stream)
        args.output.write(variant_stream.vcf_obj.get_header())
        args.output.write('\n')
        return set_missing(variant_stream, deletion_size_map, args.output, max_size, args.sr_cutoff)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
