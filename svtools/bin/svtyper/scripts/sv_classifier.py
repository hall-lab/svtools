#!/usr/bin/env python

import argparse, sys, copy, gzip, os
import math, time, re
import numpy as np
from scipy import stats
from collections import Counter
from argparse import RawTextHelpFormatter
from operator import itemgetter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.2 $"
__date__ = "$Date: 2014-04-28 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
sv_classifier.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: classify structural variants")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-g', '--gender', metavar='FILE', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
    parser.add_argument('-e', '--exclude', metavar='FILE', dest='exclude', type=argparse.FileType('r'), required=False, default=None, help='list of samples to exclude from classification algorithms')
    parser.add_argument('-a', '--annotation', metavar='BED', dest='ae_path', type=str, default=None, help='BED file of annotated elements')
    parser.add_argument('-f', '--fraction', metavar='FLOAT', dest='f_overlap', type=float, default=0.9, help='fraction of reciprocal overlap to apply annotation to variant [0.9]')
    parser.add_argument('-s', '--slope_threshold', metavar='FLOAT', dest='slope_threshold', type=float, default=1.0, help='minimum slope absolute value of regression line to classify as DEL or DUP[1.0]')
    parser.add_argument('-r', '--rsquared_threshold', metavar='FLOAT', dest='rsquared_threshold', type=float, default=0.2, help='minimum R^2 correlation value of regression line to classify as DEL or DUP [0.2]')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcf_in == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcf_in = sys.stdin
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
    def get_header(self):
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

# http://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

# test whether variant has read depth support by regression
def has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir=None):
    # slope_threshold = 0.1
    # rsquared_threshold = 0.1
    
    if 'CN' in var.active_formats:
        # allele balance list
        ab_list = []
        for s in var.sample_list:
            # if s in exclude:
            #     continue
            ab_str = var.genotype(s).get_format('AB')
            if ab_str == '.':
                ab_list.append(-1)
                continue

            ab_list.append(float(ab_str))

        # populate read-depth list, accounting for sample gender
        rd_list = []
        for s in var.sample_list:
            # if s in exclude:
            #     continue
            if (var.chrom == 'X' or var.chrom == 'Y') and gender[s] == 1:
                rd_list.append(float(var.genotype(s).get_format('CN')) * 2)
            else:
                rd_list.append(float(var.genotype(s).get_format('CN')))

        rd = np.array([ab_list, rd_list])

        # remove missing genotypes
        rd = rd[:, rd[0]!=-1]

        # ensure non-uniformity in genotype and read depth
        if len(np.unique(rd[0,:])) > 1 and len(np.unique(rd[1,:])) > 1:
            # calculate regression
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(rd)
            # print slope, intercept, r_value, var.info['SVTYPE'], var.var_id


            # write the scatterplot to a file
            if writedir is not None:
                try:
                    os.makedirs(writedir)
                except OSError as exc: # Python >2.5
                    if os.path.isdir(writedir):
                        pass
                    else: raise

                f = open('%s/reg_%s_%s_%sbp.txt' % (writedir, var.info['SVTYPE'], var.var_id, var.info['SVLEN']), 'w')
                np.savetxt(f, np.transpose(rd), delimiter='\t')
                f.close()

            if r_value ** 2 < rsquared_threshold:
                return False

            if var.info['SVTYPE'] == 'DEL':
                slope = -slope

            if slope < slope_threshold:
                return False

            return True
    return False

# test for read depth support of low frequency variants
def has_low_freq_depth_support(var, gender, exclude, writedir=None):
    mad_threshold = 2
    mad_quorum = 0.5 # this fraction of the pos. genotyped results must meet the mad_threshold
    absolute_cn_diff = 0.5
    
    hom_ref_cn = []
    het_cn = []
    hom_alt_cn = []

    for s in var.sample_list:
        if s in exclude:
            continue
        if (var.chrom == 'X' or var.chrom == 'Y') and gender[s] == 1:
            cn = float(var.genotype(s).get_format('CN')) * 2
        else:
            cn = float(var.genotype(s).get_format('CN'))

        if var.genotype(s).get_format('GT') == '0/0':
            hom_ref_cn.append(cn)
        elif var.genotype(s).get_format('GT') == '0/1':
            het_cn.append(cn)
        elif var.genotype(s).get_format('GT') == '1/1':
            hom_alt_cn.append(cn)

    if len(hom_ref_cn) > 0:
        cn_mean = np.mean(hom_ref_cn)
        cn_stdev = np.std(hom_ref_cn)
        cn_median = np.median(hom_ref_cn)
        cn_mad = mad(hom_ref_cn)
    else:
        cn_mean = None
        cn_stdev = None
        cn_median = None
        cn_mad = None

    # write the cn values to a file
    if writedir is not None:

        try:
            os.makedirs(writedir)
        except OSError as exc: # Python >2.5
            if os.path.isdir(writedir):
                pass
            else: raise

        f = open('%s/mad_%s_%s_%sbp.txt' % (writedir, var.info['SVTYPE'], var.var_id, var.info['SVLEN']), 'w')
        for cn in hom_ref_cn:
            f.write('\t'.join(map(str, [0, cn, cn_mean, cn_stdev, cn_median, cn_mad])) + '\n')
        for cn in het_cn:
            f.write('\t'.join(map(str, [1, cn, cn_mean, cn_stdev, cn_median, cn_mad])) + '\n')
        for cn in hom_alt_cn:
            f.write('\t'.join(map(str, [2, cn, cn_mean, cn_stdev, cn_median, cn_mad])) + '\n')
        f.close()

    # bail after writing out diagnostic info, if no ref samples or all ref samples
    if (len(hom_ref_cn) == 0 or
        len(het_cn + hom_alt_cn) == 0):
        return False

    # tally up the pos. genotyped samples meeting the mad_threshold
    q = 0
    for cn in het_cn + hom_alt_cn:
        resid = cn - cn_median
        if var.info['SVTYPE'] == 'DEL':
            resid = -resid
        if (resid > (cn_mad * mad_threshold) and
            resid > absolute_cn_diff):
            q += 1
    # check if meets quorum
    if float(q)/len(het_cn + hom_alt_cn) > mad_quorum:
        return True
    else:
        return False

def to_bnd(var):
    var1 = copy.deepcopy(var)
    var2 = copy.deepcopy(var)

    # update svtype
    var1.info['SVTYPE'] = 'BND'
    var2.info['SVTYPE'] = 'BND'

    # update variant id
    var1.info['EVENT'] = var.var_id
    var2.info['EVENT'] = var.var_id
    var1.var_id = var.var_id + "_1"
    var2.var_id = var.var_id + "_2"
    var1.info['MATEID'] = var2.var_id
    var2.info['MATEID'] = var1.var_id
    
    # update position
    var2.pos = var.info['END']

    # update CIPOS and CIEND
    var2.info['CIPOS'] = var.info['CIEND']
    var2.info['CIEND'] = var.info['CIPOS']
    var2.info['CIPOS95'] = var.info['CIEND95']
    var2.info['CIEND95'] = var.info['CIPOS95']

    # delete svlen and END
    del var1.info['SVLEN']
    del var2.info['SVLEN']
    del var1.info['END']
    del var2.info['END']

    # add SECONDARY to var2
    var2.info['SECONDARY'] = True

    if var.info['SVTYPE'] == 'DEL':
        var1.alt = 'N[%s:%s[' % (var.chrom, var.info['END'])
        var2.alt = ']%s:%s]N' % (var.chrom, var.pos)

    elif var.info['SVTYPE'] == 'DUP':
        var1.alt = ']%s:%s]N' % (var.chrom, var.info['END'])
        var2.alt = 'N[%s:%s[' % (var.chrom, var.pos)
    return var1, var2

def reciprocal_overlap(a, b_list):
    overlap = 0
    b_aggregate = 0
    # catch divide by zero error
    if a[1] == a[0]:
        return 0

    # update the overlap and b_aggregate
    for b in b_list:
        b_aggregate += (b[1] - b[0])
        overlap += float(min(a[1], b[1]) - max(a[0], b[0]))

    # catch divide by zero error
    if b_aggregate == 0:
        return 0

    return min(overlap / (a[1] - a[0]), overlap / b_aggregate)

def collapse_bed_records(bed_list):
    bed_list_sorted = sorted(bed_list, key=itemgetter(1))

    collapsed_bed_list = []
    
    i = 0
    curr_rec = bed_list_sorted[i]
    while i < len(bed_list_sorted):
        # end at last element in list
        if i == len(bed_list_sorted) - 1:
            collapsed_bed_list.append(copy.copy(curr_rec))
            break

        # load next entry
        next_rec = bed_list_sorted[i + 1]
        # merge is overlap
        if curr_rec[1] >= next_rec[0]:
            # print curr_rec, next_rec
            curr_rec[1] = next_rec[1]
            # print curr_rec
            i += 1
        # write out if no overlap
        else:
            collapsed_bed_list.append(copy.copy(curr_rec))
            i += 1
            curr_rec = bed_list_sorted[i]

    # print 'collapsed:', collapsed_bed_list
    return collapsed_bed_list

def annotation_intersect(var, ae_dict, threshold):
    best_frac_overlap = 0
    best_feature = ''
    slop = 0

    # dictionary with number of bases of overlap for each class
    class_overlap = {}
    
    # first check for reciprocal overlap
    if var.chrom in ae_dict:
        var_start = var.pos
        var_end = int(var.info['END'])
        i = 0
        while 1:
            # bail if end of dict
            if i >= len(ae_dict[var.chrom]):
                break
            feature = ae_dict[var.chrom][i]
            if feature[0] - slop < var_end:
                if feature[1] + slop > var_start:
                    try:
                        class_overlap[feature[2]].append(feature)
                    except KeyError:
                        class_overlap[feature[2]] = [feature]
            else:
                break
            i += 1

        # print class_overlap
        for me_class in class_overlap:
            class_overlap[me_class] = collapse_bed_records(class_overlap[me_class])
            # print 'class_overlap[me_class]:', class_overlap[me_class]
            # print 'recip:', reciprocal_overlap([var_start, var_end], class_overlap[me_class])

            frac_overlap = reciprocal_overlap([var_start, var_end], class_overlap[me_class])
            if frac_overlap > best_frac_overlap:
                best_frac_overlap = frac_overlap
                best_feature = me_class


        if best_frac_overlap >= threshold:
            return best_feature

    return None

# primary function
def sv_classify(vcf_in, gender_file, exclude_file, ae_dict, f_overlap, slope_threshold, rsquared_threshold):
    vcf_out = sys.stdout
    vcf = Vcf()
    header = []
    in_header = True
    min_pos_samps_for_regression = 10

    gender = {}
    # read sample genders
    for line in gender_file:
        v = line.rstrip().split('\t')
        gender[v[0]] = int(v[1])

    exclude = []
    if exclude_file is not None:
        for line in exclude_file:
            exclude.append(line.rstrip())

    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                continue
            else:
                in_header = False
                vcf.add_header(header)
                # write the output header
                vcf_out.write(vcf.get_header() + '\n')

        # split variant line, quick pre-check if the SVTYPE is BND, and skip if so
        v = line.rstrip().split('\t')

        info = v[7].split(';')
        svtype = None
        for x in info:
            if x.startswith('SVTYPE='):
                svtype = x.split('=')[1]
                break

        # bail if not DEL or DUP prior to reclassification
        if svtype not in ['DEL', 'DUP']:
            vcf_out.write(line)
            continue

        # parse the VCF line
        var = Variant(v, vcf)

        # check intersection with mobile elements
        if ae_dict is not None and var.info['SVTYPE'] in ['DEL']:
            ae = annotation_intersect(var, ae_dict, f_overlap)
            if ae is not None:
                if ae.startswith('SINE') or ae.startswith('LINE') or ae.split('|')[2].startswith('SVA'):
                    ae = 'ME:' + ae
                var.alt = '<DEL:%s>' % ae
                var.info['SVTYPE'] = 'MEI'
                vcf_out.write(var.get_var_string() + '\n')
                continue

        # # write to directory
        # writedir = 'data/r11.100kb.dup'

        # annotate based on read depth
        if var.info['SVTYPE'] in ['DEL', 'DUP']:
            # count the number of positively genotyped samples
            num_pos_samps = 0;
            for s in var.sample_list:
                if s in exclude:
                    continue
                if var.genotype(s).get_format('GT') not in ["./.", "0/0"]:
                    num_pos_samps += 1

            if num_pos_samps < min_pos_samps_for_regression:
                if has_low_freq_depth_support(var, gender, exclude):
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/low_freq_rd')
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/low_freq_rd')
                    # write variant
                    vcf_out.write(var.get_var_string() + '\n')
                else:
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/low_freq_no_rd')
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/low_freq_no_rd')
                    for m_var in to_bnd(var):
                        vcf_out.write(m_var.get_var_string() + '\n')
            else:
                if has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold):
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/high_freq_rd')
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/high_freq_rd')
                    # write variant
                    vcf_out.write(var.get_var_string() + '\n')
                else:
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/high_freq_no_rd')
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/high_freq_no_rd')
                    for m_var in to_bnd(var):
                        vcf_out.write(m_var.get_var_string() + '\n')
    vcf_out.close()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # load the annotated elements
    ae_dict = None
    if args.ae_path is not None:
        if args.ae_path.endswith('.gz'):
            ae_bedfile = gzip.open(args.ae_path, 'rb')
        else:
            ae_bedfile = open(args.ae_path, 'r')
        ae_dict = {}

        for line in ae_bedfile:
            v = line.rstrip().split('\t')
            if len(v) < 4:
                continue
            v[1] = int(v[1])
            v[2] = int(v[2])
            if v[0] in ae_dict:
                ae_dict[v[0]].append(v[1:])
            else:
                ae_dict[v[0]] = [v[1:]]

    # call primary function
    sv_classify(args.vcf_in,
                args.gender,
                args.exclude,
                ae_dict,
                args.f_overlap,
                args.slope_threshold,
                args.rsquared_threshold)

    # close the files
    args.vcf_in.close()
    args.gender.close()
    args.exclude.close()
    if args.ae_path is not None:
        ae_bedfile.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
