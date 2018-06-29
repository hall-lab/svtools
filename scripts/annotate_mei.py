#!/usr/bin/env python

import argparse
import sys

import svtools.utils as su
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from svtools.sv_classifier import annotation_intersect, get_ae_dict

# NOTE This is duplicated from filter_del and should be moved into the main repo and used to deduplicate a bunch of nonsense

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Standalone script to annotate MEIs independent of reclassification")
    parser.add_argument('-i', '--input', metavar='<VCF>', default=None, help='VCF input')
    parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
    parser.add_argument('-a', '--annotation', metavar='<BED>', dest='ae_path', type=str, required=True, help='BED file of annotated elements')
    parser.add_argument('-f', '--fraction', metavar='<FLOAT>', dest='f_overlap', type=float, default=0.9, help='fraction of reciprocal overlap to apply annotation to variant [0.9]')
    args = parser.parse_args()

    ae_dict = get_ae_dict(args.ae_path)

    with su.InputStream(args.input) as stream:
        variant_stream = VCFReader(stream)
        args.vcf_out.write(variant_stream.vcf_obj.get_header())
        args.vcf_out.write('\n')
        for v in variant_stream:
            svtype = v.get_info('SVTYPE')
            if svtype == 'DEL':
                ae = annotation_intersect(v, ae_dict, args.f_overlap)
                if ae is not None:
                    if ae.startswith('SINE') or ae.startswith('LINE') or ae.split('|')[2].startswith('SVA'):
                        ae = 'ME:' + ae
                    v.alt = '<DEL:{}>'.format(ae)
                    v.set_info('SVTYPE', 'MEI')
            args.vcf_out.write(v.get_var_string(use_cached_gt_string=True))
            args.vcf_out.write('\n')
