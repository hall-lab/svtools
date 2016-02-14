import svtools.l_bp as l_bp

import sys
import os
import heapq
import argparse
from tempfile import gettempdir
from collections import namedtuple

Keyed = namedtuple("Keyed", ["key", "obj"])
def merge(*iterables):
   keyed_iterables = [(Keyed(l_bp.vcf_line_key(obj), obj) for obj in iterable)
                        for iterable in iterables]
   for element in heapq.merge(*keyed_iterables):
       yield element.obj

class Lsort(object):
    def __init__(self, vcf_file_names, tempdir=None, batchsize=200):
        if tempdir:
            self.tempdir = tempdir
        else:
            self.tempdir = gettempdir()
        self.batchsize = batchsize
        self.vcf_file_names = vcf_file_names
        self.vcf_lines = []
        self.vcf_headers = []
        self.temp_files = []

    def execute(self):
        
        counter = 0
        for vcf_file_name in self.vcf_file_names:
            samples = l_bp.parse_vcf(vcf_file_name, self.vcf_lines, self.vcf_headers)
            for sample in samples:
                self.vcf_headers.append("##SAMPLE=<ID=" + sample + ">\n")
            counter += 1
            if counter > self.batchsize:
                self.vcf_lines.sort(key=l_bp.vcf_line_key)
                self.write_temp_file()
                counter = 0
        # no need to write the final batch to file
        self.write_header()

        self.vcf_lines.sort(key=l_bp.vcf_line_key)
        iterables = self.temp_files + [self.vcf_lines]
        sys.stdout.writelines(merge(*iterables))

    def close_tempfiles(self):
        for tmp in self.temp_files:
            tmp.close()
            os.remove(tmp.name)

    def write_header(self):
        self.vcf_headers.append("##INFO=<ID=SNAME,Number=.,Type=String," + \
            "Description=\"Source sample name\">\n")
        self.vcf_headers.append("##INFO=<ID=ALG,Number=1,Type=String," + \
            "Description=\"Evidence PDF aggregation algorithm\">\n")
        self.vcf_headers.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" + \
            "VARIOUS\n")
        self.vcf_headers.sort(cmp=l_bp.header_line_cmp)
        sys.stdout.writelines(self.vcf_headers)

    def write_temp_file(self):
        temp_outfile = open(os.path.join(self.tempdir,'%06i'%len(self.temp_files)),'w+b',64*1024)
        temp_outfile.writelines(self.vcf_lines)
        temp_outfile.flush()
        temp_outfile.seek(0)
        self.temp_files.append(temp_outfile)

        #vcf_line array
        self.vcf_lines = []

def description():
    return 'sort N LUMPY VCF files into a single file'

def epilog():
    return '''Specify -t to override where temporary files are placed. Use -b to control the amount of memory required. 
    This will vary depending on the number of lines in your input files.'''

def add_arguments_to_parser(parser):
    parser.add_argument('vcf_files', metavar='<VCF>', nargs='+', help='VCF files to combine and sort')
    parser.add_argument('-t', '--tempdir', metavar='<DIRECTORY_PATH>', default=gettempdir(), help='temporary directory')
    parser.add_argument('-b', '--batchsize', metavar='<INT>', type=int, default=200, help='number of files to sort in batch')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description(), epilog=epilog())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    sorter = Lsort(args.vcf_files, tempdir=args.tempdir, batchsize=args.batchsize)
    sorter.execute()

if __name__ == "__main__":
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
