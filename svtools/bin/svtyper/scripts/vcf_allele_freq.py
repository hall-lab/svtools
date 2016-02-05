#!/usr/bin/env python

import pysam
import argparse, sys
import math, time, re
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-04-22 09:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
vcf_allele_freq.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Add allele frequency information to a VCF file")
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
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]

    # return the VCF header
    def get_header(self, include_samples=True):
        if include_samples:
            header = '\n'.join(['##fileformat=' + self.file_format,
                                '##fileDate=' + time.strftime('%Y%m%d'),
                                '##reference=' + self.reference] + \
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

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

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
        if len(self.active_formats) == 0:
            s = '\t'.join(map(str,[
                self.chrom,
                self.pos,
                self.var_id,
                self.ref,
                self.alt,
                '%0.2f' % self.qual,
                self.filter,
                self.get_info_string()
            ]))
        else:
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
def add_af(vcf_file):
    in_header = True
    header = []
    breakend_dict = {} # cache to hold unmatched generic breakends for genotyping
    vcf = Vcf()
    vcf_out = sys.stdout

    # read input VCF
    for line in vcf_file:
        if in_header:
            if line.startswith('##'):
                header.append(line) 
                continue
            elif line.startswith('#CHROM'):
                v = line.rstrip().split('\t')
                header.append('\t'.join(v[:8]))

                in_header = False
                vcf.add_header(header)
                
                vcf.add_info('AF', 'A', 'Float', 'Allele Frequency, for each ALT allele, in the same order as listed')
                vcf.add_info('NSAMP', '1', 'Integer', 'Number of samples with non-reference genotypes')

                # write header
                vcf_out.write(vcf.get_header(include_samples=False))
                vcf_out.write('\t' + '\t'.join(v[8:]) + '\n')
            continue

        v = line.rstrip().split('\t')
        var = Variant(v[:8], vcf)

        # extract genotypes from VCF
        num_alt = len(var.alt.split(','))
        alleles = [0] * (num_alt + 1)
        num_samp = 0

        for i in xrange(9,len(v)):
            gt_string = v[i].split(':')[0]

            if '.' in  gt_string:
                continue
            gt = gt_string.split('/')
            if len(gt) == 1:
                gt = gt_string.split('|')
            gt = map(int, gt)

            for i in xrange(len(gt)):
                alleles[gt[i]] += 1

            # iterate the number of non-reference samples
            if sum(gt) > 0:
                num_samp += 1

        allele_sum = float(sum(alleles))
        allele_freq = ['.'] * len(alleles)

        # populate AF
        if allele_sum > 0:
            for i in xrange(len(alleles)):
                allele_freq[i] = alleles[i] / allele_sum
            var.info['AF'] = ','.join(map(str, ['%.4g' % a for a in allele_freq[1:]]))
        else:
            var.info['AF'] = ','.join(map(str, allele_freq[1:]))
        
        # populate NSAMP
        var.info['NSAMP'] = num_samp

        # after all samples have been processed, write
        vcf_out.write(var.get_var_string()
                      + '\t'
                      + '\t'.join(v[8:])
                      + '\n')
    vcf_out.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    add_af(args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
