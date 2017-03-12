#!/usr/bin/env python

import pysam
import argparse, sys
import math, time, re
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-09-28 11:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
vcf_modify_header.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Add or remove lines from the header")
    parser.add_argument('-i', '--id', dest='vcf_id', metavar='STR', type=str, required=False, help='field id')
    parser.add_argument('-c', '--category', dest='category', metavar='STR', type=str, help='INFO, FORMAT, FILTER')
    parser.add_argument('-t', '--type', dest='type', metavar='STR', type=str, required=False, help='Number, String, Float, Integer')
    parser.add_argument('-n', '--number', dest='number', metavar='STR', type=str, required=False, help='integer, A, R, G, or .')
    parser.add_argument('-d', '--description', dest='description', metavar='STR', type=str, required=False, help='description')
    parser.add_argument(metavar='vcf', dest='input_vcf', nargs='?', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input_vcf = sys.stdin
    # send back the user input
    return args

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        # self.fasta = fasta
        self.reference = ''
        self.sample_list = []
        self.info_list = []
        self.format_list = []
        self.filter_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')

    def add_header(self, header):
        for line in header:
            if line.split('=')[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##INFO':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_info(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##ALT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_alt(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FORMAT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_format(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FILTER':
                a = line[line.find('<')+1:-2]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_filter(*[b.split('=')[1] for b in r.findall(a)])
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]

    # return the VCF header
    def get_header(self, include_samples=True):
        if include_samples:
            header = '\n'.join(['##fileformat=' + self.file_format,
                                '##fileDate=' + time.strftime('%Y%m%d'),
                                '##reference=' + self.reference] + \
                               [f.hstring for f in self.filter_list] + \
                               [i.hstring for i in self.info_list] + \
                               [a.hstring for a in self.alt_list] + \
                               [f.hstring for f in self.format_list] + \
                               ['\t'.join([
                                   '#CHROM',
                                   'POS',
                                   'ID',
                                   'REF',
                                   'ALT',
                                   'QUAL',
                                   'FILTER',
                                   'INFO',
                                   'FORMAT'] + \
                                          self.sample_list
                                      )])
        else:
            header = '\n'.join(['##fileformat=' + self.file_format,
                                '##fileDate=' + time.strftime('%Y%m%d'),
                                '##reference=' + self.reference] + \
                               [f.hstring for f in self.filter_list] + \
                               [i.hstring for i in self.info_list] + \
                               [a.hstring for a in self.alt_list] + \
                               [f.hstring for f in self.format_list] + \
                               ['\t'.join([
                                   '#CHROM',
                                   'POS',
                                   'ID',
                                   'REF',
                                   'ALT',
                                   'QUAL',
                                   'FILTER',
                                   'INFO']
                                          )])
        return header

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def remove_info(self, id):
        for i in self.info_list:
            if i.id == id:
                self.info_list.remove(i)
                return

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_filter(self, id, desc):
        if id not in [f.id for f in self.filter_list]:
            filter = self.Filter(id, desc)
            self.filter_list.append(filter)

    def add_sample(self, name):
        self.sample_list.append(name)

    # get the VCF column index of a sample
    # NOTE: this is zero-based, like python arrays
    def sample_to_col(self, sample):
        return self.sample_list.index(sample) + 9

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

    class Filter(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FILTER=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

class Variant(object):
    def __init__(self, var_list, vcf):
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        if var_list[5] == '.':
            self.qual = 0
        else:
            self.qual = float(var_list[5])
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.active_formats = list()
        self.gts = dict()
        
        # fill in empty sample genotypes
        if len(var_list) < 8:
            sys.stderr.write('\nError: VCF file must have at least 8 columns\n')
            exit(1)
        if len(var_list) < 9:
            var_list.append("GT")

        # make a genotype for each sample at variant
        for s in self.sample_list:
            try:
                s_gt = var_list[vcf.sample_to_col(s)].split(':')[0]
                self.gts[s] = Genotype(self, s, s_gt)
                # import the existing fmt fields
                for j in zip(var_list[8].split(':'), var_list[vcf.sample_to_col(s)].split(':')):
                    self.gts[s].set_format(j[0], j[1])
            except IndexError:
                self.gts[s] = Genotype(self, s, './.')

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        return self.info[field]

    def get_info_string(self):
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)

    def get_format_string(self):
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)

    def genotype(self, sample_name):
        if sample_name in self.sample_list:
            return self.gts[sample_name]
        else:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')

    def get_var_string(self):
        s = '\t'.join(map(str,[
            self.chrom,
            self.pos,
            self.var_id,
            self.ref,
            self.alt,
            '%0.2f' % self.qual,
            self.filter,
            self.get_info_string(),
            self.get_format_string(),
            '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        ]))
        return s

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.append(field)
                # sort it to be in the same order as the format_list in header
                self.variant.active_formats.sort(key=lambda x: [f.id for f in self.variant.format_list].index(x))
        else:
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
        g_list = list()
        for f in self.variant.active_formats:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(self.format[f])
            else:
                g_list.append('.')
        return ':'.join(map(str,g_list))

# primary function
def mod_header(vcf_id,
               category,
               f_type,
               number,
               description,
               vcf_file):
    in_header = True
    header = []
    breakend_dict = {} # cache to hold unmatched generic breakends for genotyping
    vcf = Vcf()
    vcf_out = sys.stdout

    # read input VCF
    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line) 
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)
                if category == 'INFO':
                    vcf.remove_info(vcf_id)
                    vcf.add_info(vcf_id, number, f_type, description)
                elif category == 'FORMAT':
                    vcf.add_format(vcf_id, number, f_type, description)
                elif category == 'FILTER':
                    vcf.add_filter(vcf_id, description)

                # write the output header
                if len(vcf_samples) > 0:
                    vcf_out.write(vcf.get_header(include_samples=True) + '\n')
                else:
                    vcf_out.write(vcf.get_header(include_samples=False) + '\n')
        vcf_out.write(line)
    vcf_out.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    mod_header(args.vcf_id,
               args.category,
               args.type,
               args.number,
               args.description,
               args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
