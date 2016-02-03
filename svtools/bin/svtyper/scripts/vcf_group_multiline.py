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
vcf_group_multiline.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group multiline variants prior vcf_paste.py")
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
def sv_genotype(vcf_file):
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
                # if detailed:
                vcf.add_format('GQ', 1, 'Float', 'Genotype quality')
                vcf.add_format('SQ', 1, 'Float', 'Phred-scaled probability that this site is variant (non-reference in this sample')
                vcf.add_format('GL', 'G', 'Float', 'Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible gen\
otype generated from the reference and alternate alleles given the sample ploidy')
                vcf.add_format('DP', 1, 'Integer', 'Read depth')
                vcf.add_format('RO', 1, 'Integer', 'Reference allele observation count, with partial observations recorded fractionally')
                vcf.add_format('AO', 'A', 'Integer', 'Alternate allele observations, with partial observations recorded fractionally')
                vcf.add_format('QR', 1, 'Integer', 'Sum of quality of reference observations')
                vcf.add_format('QA', 'A', 'Integer', 'Sum of quality of alternate observations')
                vcf.add_format('RS', 1, 'Integer', 'Reference allele split-read observation count, with partial observations recorded fractionally')
                vcf.add_format('AS', 'A', 'Integer', 'Alternate allele split-read observation count, with partial observations recorded fractionally')
                vcf.add_format('RP', 1, 'Integer', 'Reference allele paired-end observation count, with partial observations recorded fractionally')
                vcf.add_format('AP', 'A', 'Integer', 'Alternate allele paired-end observation count, with partial observations recorded fractionally')
                vcf.add_format('AB', 'A', 'Float', 'Allele balance, fraction of observations from alternate allele, QA/(QR+QA)')

                # write the output header
                if len(vcf_samples) > 0:
                    vcf_out.write(vcf.get_header(include_samples=True) + '\n')
                else:
                    vcf_out.write(vcf.get_header(include_samples=False) + '\n')

        v = line.rstrip().split('\t')
        var = Variant(v, vcf)

        # genotype generic breakends
        if var.info['SVTYPE']=='BND':
            if var.info['MATEID'] in breakend_dict:
                var2 = var
                var = breakend_dict[var.info['MATEID']]
                chromA = var.chrom
                chromB = var2.chrom
                posA = var.pos
                posB = var2.pos
                # confidence intervals
                ciA = [posA + ci for ci in map(int, var.info['CIPOS'].split(','))]
                ciB = [posB + ci for ci in map(int, var2.info['CIPOS'].split(','))]

                # infer the strands from the alt allele
                if var.alt[-1] == '[' or var.alt[-1] == ']':
                    o1 = '+'
                else: o1 = '-'
                if var2.alt[-1] == '[' or var2.alt[-1] == ']':
                    o2 = '+'
                else: o2 = '-'
            else:
                breakend_dict[var.var_id] = var
                continue
        else:
            chromA = var.chrom
            chromB = var.chrom
            posA = var.pos
            posB = int(var.get_info('END'))
            # confidence intervals
            ciA = [posA + ci for ci in map(int, var.info['CIPOS'].split(','))]
            ciB = [posB + ci for ci in map(int, var.info['CIEND'].split(','))]
            if var.get_info('SVTYPE') == 'DEL':
                o1, o2 =  '+', '-'
            elif var.get_info('SVTYPE') == 'DUP':
                o1, o2 =  '-', '+'
            elif var.get_info('SVTYPE') == 'INV':
                o1, o2 =  '+', '+'

        # # increment the negative strand values (note position in VCF should be the base immediately left of the breakpoint junction)
        # if o1 == '-': posA += 1
        # if o2 == '-': posB += 1
        # # if debug: print posA, posB

        # # for i in xrange(len(bam_list)):
        # for sample in sample_list:
        #     '''
        #     Breakend A
        #     '''
        #     # Count splitters
        #     ref_counter_a = Counter()
        #     spl_counter_a = Counter()
        #     ref_scaled_counter_a = Counter()
        #     spl_scaled_counter_a = Counter()

        #     for ref_read in sample.bam.fetch(chromA, max(posA - padding, 0), posA + padding + 1):
        #         if not ref_read.is_duplicate and not ref_read.is_unmapped:
        #             for p in xrange(ref_read.pos + 1, ref_read.aend + 1):
        #                 if p - ref_read.pos >= splflank and ref_read.aend - p >= splflank:
        #                     ref_counter_a[p] += 1
        #                     ref_scaled_counter_a[p] += (1-10**(-ref_read.mapq/10.0))
        #     for spl_read in sample.spl_bam.fetch(chromA, max(posA - padding, 0), posA + padding + 1):
        #         if not spl_read.is_duplicate and not spl_read.is_unmapped:
        #             if o1 == '+' and spl_read.cigar[0][0] == 0:
        #                 # if debug: print 'o1+', spl_read.aend
        #                 spl_counter_a[spl_read.aend] += 1
        #                 spl_scaled_counter_a[spl_read.aend] += (1-10**(-spl_read.mapq/10.0))
        #             elif o1 == '-' and spl_read.cigar[-1][0] == 0:
        #                 # if debug: print 'o1-', spl_read.pos + 1
        #                 spl_counter_a[spl_read.pos + 1] += 1
        #                 spl_scaled_counter_a[spl_read.pos + 1] += (1-10**(-spl_read.mapq/10.0))

        #     # Count paired-end discordant and concordants
        #     (conc_counter_a,
        #      disc_counter_a,
        #      conc_scaled_counter_a,
        #      disc_scaled_counter_a) = count_pairedend(chromA, posA, ciA,
        #                                               chromB, posB, ciB,
        #                                               o1, o2,
        #                                               var.info['SVTYPE'],
        #                                               sample,
        #                                               z, discflank)
        #     '''
        #     Breakend B
        #     '''
        #     # Count splitters
        #     ref_counter_b = Counter()
        #     spl_counter_b = Counter()
        #     ref_scaled_counter_b = Counter()
        #     spl_scaled_counter_b = Counter()

        #     for ref_read in sample.bam.fetch(chromB, max(posB - padding, 0), posB + padding + 1):
        #         if not ref_read.is_duplicate and not ref_read.is_unmapped:
        #             for p in xrange(ref_read.pos + 1, ref_read.aend + 1):
        #                 if p - ref_read.pos >= splflank and ref_read.aend - p >= splflank:
        #                     ref_counter_b[p] += 1
        #                     ref_scaled_counter_b[p] += (1-10**(-ref_read.mapq/10.0))
        #     for spl_read in sample.spl_bam.fetch(chromB, max(posB - padding, 0), posB + padding + 1):
        #         if not spl_read.is_duplicate and not spl_read.is_unmapped:
        #             if o2 == '+' and spl_read.cigar[0][0] == 0:
        #                 spl_counter_b[spl_read.aend] += 1
        #                 # if debug: print 'o2+', spl_read.aend
        #                 spl_scaled_counter_b[spl_read.aend] += (1-10**(-spl_read.mapq/10.0))
        #             elif o2 == '-' and spl_read.cigar[-1][0] == 0:
        #                 # if debug: print 'o2-', spl_read.pos + 1
        #                 spl_counter_b[spl_read.pos + 1] += 1
        #                 spl_scaled_counter_b[spl_read.pos + 1] += (1-10**(-spl_read.mapq/10.0))
            
        #     # tally up the splitters
        #     sr_ref_a = int(round(sum(ref_counter_a[p] for p in xrange(posA - split_slop, posA + split_slop + 1)) / float(2 * split_slop + 1)))
        #     sr_spl_a = sum(spl_counter_a[p] for p in xrange(posA-split_slop, posA+split_slop + 1))
        #     sr_ref_b = int(round(sum(ref_counter_b[p] for p in xrange(posB - split_slop, posB + split_slop + 1)) / float(2 * split_slop + 1)))
        #     sr_spl_b = sum(spl_counter_b[p] for p in xrange(posB - split_slop, posB + split_slop + 1))

        #     sr_ref_scaled_a = sum(ref_scaled_counter_a[p] for p in xrange(posA - split_slop, posA + split_slop + 1)) / float(2 * split_slop + 1)
        #     sr_spl_scaled_a = sum(spl_scaled_counter_a[p] for p in xrange(posA-split_slop, posA+split_slop + 1))
        #     sr_ref_scaled_b = sum(ref_scaled_counter_b[p] for p in xrange(posB - split_slop, posB + split_slop + 1)) / float(2 * split_slop + 1)
        #     sr_spl_scaled_b = sum(spl_scaled_counter_b[p] for p in xrange(posB - split_slop, posB + split_slop + 1))

        #     # Count paired-end discordants and concordants
        #     (conc_counter_b,
        #      disc_counter_b,
        #      conc_scaled_counter_b,
        #      disc_scaled_counter_b) = count_pairedend(chromB, posB, ciB,
        #                                               chromA, posA, ciA,
        #                                               o2, o1,
        #                                               var.info['SVTYPE'],
        #                                               sample,
        #                                               z, discflank)
        #     if debug:
        #         print '--------------------'
        #         print sample.name
        #         print 'sr_a', '(ref, alt)', sr_ref_a, sr_spl_a
        #         print 'pe_a', '(ref, alt)', conc_counter_a, disc_counter_a
        #         print 'sr_b', '(ref, alt)', sr_ref_b, sr_spl_b
        #         print 'pe_b', '(ref, alt)', conc_counter_b, disc_counter_b
        #         print 'sr_a_scaled', '(ref, alt)', sr_ref_scaled_a, sr_spl_scaled_a
        #         print 'pe_a_scaled', '(ref, alt)', conc_scaled_counter_a, disc_scaled_counter_a
        #         print 'sr_b_scaled', '(ref, alt)', sr_ref_scaled_b, sr_spl_scaled_b
        #         print 'pe_b_scaled', '(ref, alt)', conc_scaled_counter_b, disc_scaled_counter_b

        #     # merge the breakend support
        #     split_ref = 0 # set these to zero unless there are informative alt bases for the ev type
        #     disc_ref = 0
        #     split_alt = sr_spl_a + sr_spl_b
        #     if split_alt > 0:
        #         split_ref = sr_ref_a + sr_ref_b
        #     disc_alt = disc_counter_a + disc_counter_b
        #     if disc_alt > 0:
        #         disc_ref = conc_counter_a + conc_counter_b
        #     if split_alt == 0 and disc_alt == 0:
        #         split_ref = sr_ref_a + sr_ref_b
        #         disc_ref = conc_counter_a + conc_counter_b

        #     split_scaled_ref = 0 # set these to zero unless there are informative alt bases for the ev type
        #     disc_scaled_ref = 0
        #     split_scaled_alt = sr_spl_scaled_a + sr_spl_scaled_b
        #     if int(split_scaled_alt) > 0:
        #         split_scaled_ref = sr_ref_scaled_a + sr_ref_scaled_b
        #     disc_scaled_alt = disc_scaled_counter_a + disc_scaled_counter_b
        #     if int(disc_scaled_alt) > 0:
        #         disc_scaled_ref = conc_scaled_counter_a + conc_scaled_counter_b
        #     if int(split_scaled_alt) == 0 and int(disc_scaled_alt) == 0: # if no alt alleles, set reference
        #         split_scaled_ref = sr_ref_scaled_a + sr_ref_scaled_b
        #         disc_scaled_ref = conc_scaled_counter_a + conc_scaled_counter_b

        #     if split_scaled_alt + split_scaled_ref + disc_scaled_alt + disc_scaled_ref > 0:
        #         # get bayesian classifier
        #         if var.info['SVTYPE'] == "DUP": is_dup = True
        #         else: is_dup = False
        #         gt_lplist = bayes_gt(int(split_weight * split_scaled_ref) + int(disc_weight * disc_scaled_ref), int(split_weight * split_scaled_alt) + int(disc_weight * disc_scaled_alt), is_dup)
        #         gt_idx = gt_lplist.index(max(gt_lplist))

        #         # print log probabilities of homref, het, homalt
        #         if debug:
        #             print gt_lplist

        #         # set the overall variant QUAL score and sample specific fields
        #         var.genotype(sample.name).set_format('GL', ','.join(['%.0f' % x for x in gt_lplist]))
        #         var.genotype(sample.name).set_format('DP', int(split_scaled_ref + split_scaled_alt + disc_scaled_ref + disc_scaled_alt))
        #         var.genotype(sample.name).set_format('AO', int(split_scaled_alt + disc_scaled_alt))
        #         var.genotype(sample.name).set_format('RO', int(split_scaled_ref + disc_scaled_ref))
        #         # if detailed:
        #         var.genotype(sample.name).set_format('AS', int(split_scaled_alt))
        #         var.genotype(sample.name).set_format('RS', int(split_scaled_ref))
        #         var.genotype(sample.name).set_format('AP', int(disc_scaled_alt))
        #         var.genotype(sample.name).set_format('RP', int(disc_scaled_ref))

        #         # assign genotypes
        #         gt_sum = 0
        #         for gt in gt_lplist:
        #             try:
        #                 gt_sum += 10**gt
        #             except OverflowError:
        #                 gt_sum += 0
        #         if gt_sum > 0:
        #             gt_sum_log = math.log(gt_sum, 10)
        #             sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
        #             if 1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log) == 0:
        #                 phred_gq = 200                    
        #             else:
        #                 phred_gq = abs(-10 * math.log(1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log), 10))
        #             var.genotype(sample.name).set_format('GQ', phred_gq)
        #             var.genotype(sample.name).set_format('SQ', sample_qual)
        #             var.qual += sample_qual
        #             if gt_idx == 1:
        #                 var.genotype(sample.name).set_format('GT', '0/1')
        #             elif gt_idx == 2:
        #                 var.genotype(sample.name).set_format('GT', '1/1')
        #             elif gt_idx == 0:
        #                 var.genotype(sample.name).set_format('GT', '0/0')
        #         else:
        #             var.genotype(sample.name).set_format('GQ', '.')
        #             var.genotype(sample.name).set_format('SQ', '.')
        #             var.genotype(sample.name).set_format('GT', './.')
        #     else:
            # var.genotype(sample.name).set_format('GT', './.')
            # var.qual = 0
            # var.genotype(sample.name).set_format('GQ', '.')
            # var.genotype(sample.name).set_format('GL', '.')
            # var.genotype(sample.name).set_format('DP', 0)
            # var.genotype(sample.name).set_format('AO', 0)
            # var.genotype(sample.name).set_format('RO', 0)
            # # if detailed:
            # var.genotype(sample.name).set_format('AS', 0)
            # var.genotype(sample.name).set_format('RS', 0)
            # var.genotype(sample.name).set_format('AP', 0)
            # var.genotype(sample.name).set_format('RP', 0)

        # after all samples have been processed, write
        vcf_out.write(var.get_var_string() + '\n')
        if var.info['SVTYPE'] == 'BND':
            var2.qual = var.qual
            var2.active_formats = var.active_formats
            var2.genotype = var.genotype
            vcf_out.write(var2.get_var_string() + '\n')
    vcf_out.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    sv_genotype(args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
