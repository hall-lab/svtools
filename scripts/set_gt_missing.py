#!/usr/bin/env python

import argparse, sys, copy, gzip, time, math, re
from collections import Counter, defaultdict, namedtuple
import statsmodels.formula.api as smf
from operator import itemgetter
import warnings
from svtools.vcf.file import Vcf
from svtools.vcf.genotype import Genotype
from svtools.vcf.variant import Variant
import svtools.utils as su

inssz_rec = namedtuple ('inssz_rec', 'inssz_mean inssz_sd')

def set_gt_missing(vcf_in, vcf_out, inssz_file, maxlen):

    vcf = Vcf()
    header = []
    in_header = True
    inssz={}
    
    for line in inssz_file:
        v = line.rstrip().split('\t')
        inssz[v[0]] = inssz_rec(float(v[1]), float(v[2]))
    
    for line in vcf_in:
        if in_header:
            if line[0] == "#":
               header.append(line)
               continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf.add_info('GTMOD', 1, 'String', 'Small del, may be undetectable in some samples')
                vcf.add_format('GT0', 1, 'String', 'Original genotype')
                vcf_out.write(vcf.get_header() + '\n')

        v = line.rstrip().split('\t')
        info = v[7].split(';')

        var = Variant(v, vcf)
        #sys.stderr.write(var.get_info('SVTYPE')+"\t"+var.get_info('SVLEN')+"\t"+var.get_info('SR')+"\n")
        if var.get_info('SVTYPE') not in ['DEL'] or abs(int(var.get_info('SVLEN'))) > maxlen or int(var.get_info('SR'))>0:
            #sys.stderr.write("no change\n")
            vcf_out.write(line)
            continue
        
        del_len=abs(int(var.get_info('SVLEN')))
    
        changed = False
        for s in var.sample_list:
          minsz=4*inssz[s].inssz_sd
          if minsz>del_len:
              var.genotype(s).set_format("GT0", var.genotype(s).get_format("GT"))
              var.genotype(s).set_format("GT", './.')
              changed = True

        if changed:
            var.set_info('GTMOD', "True")

        vcf_out.write(var.get_var_string(use_cached_gt_string=False) + '\n')

    vcf_out.close()
    vcf_in.close()
    inssz_file.close()
    return

            



def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
    parser.add_argument('-l', '--maxlen', metavar='<INT>', dest='maxlen', type=int, default=300, required=False, help='max deletion length to consider for missingness')
    parser.add_argument('-s', '--inssz_file', metavar='<STRING>', dest='inssz_file', type=argparse.FileType('r'), default=None, required=False, help='tab-delimited insert size file, col1=sample, col2=mean inssz, col3=sd inssz')
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'refine genotypes by clustering'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.vcf_in) as stream:
        set_gt_missing(stream, args.vcf_out, args.inssz_file, args.maxlen)


if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
