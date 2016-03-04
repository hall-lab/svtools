#!/usr/bin/env python

import argparse, sys, copy, gzip, os
import math, time, re
import numpy as np
from scipy import stats
from collections import Counter
from argparse import RawTextHelpFormatter
from operator import itemgetter
from svtools.vcf.file import Vcf
from svtools.vcf.genotype import Genotype
from svtools.vcf.variant import Variant
import svtools.utils as su


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
    args.vcf_in = su.InputStream(args.vcf_in)
    # send back the user input
    return args


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

def to_bnd_strings(var):

    old_type = var.info['SVTYPE']
    old_id = var.var_id
    old_pos = var.pos
    old_end = var.info['END']
    old_ciend = var.info['CIEND']
    old_cipos = var.info['CIPOS']
    old_cipos95 = var.info['CIPOS95']
    old_ciend95 = var.info['CIEND95']


    #for both ends
    var.info['SVTYPE'] = 'BND'
    var.info['EVENT'] = old_id
    del var.info['SVLEN']
    del var.info['END']

    #var1
    var.var_id = old_id + "_1"
    var.info['MATEID'] = old_id + "_2"
    if old_type == 'DEL':
        var.alt = 'N[%s:%s[' % (var.chrom, old_end)
    else:
        var.alt = ']%s:%s]N' % (var.chrom, old_end)
    var1=var.get_var_string(True)

    #var2
    var.var_id = old_id + "_2"
    var.info['MATEID'] = old_id + "_1"
    var.info['CIPOS'] = old_ciend
    var.info['CIEND'] = old_cipos
    var.info['CIPOS95'] = old_ciend95
    var.info['CIEND95'] = old_cipos95
    var.pos = old_end
    var.info['SECONDARY'] = True
    if old_type == 'DEL':
        var.alt = ']%s:%s]N' % (var.chrom, old_pos)
    else:
        var.alt = 'N[%s:%s[' % (var.chrom, old_pos)
    var2=var.get_var_string(True)
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
        var = Variant(v, vcf, True)

        # check intersection with mobile elements
        if ae_dict is not None and var.info['SVTYPE'] in ['DEL']:
            ae = annotation_intersect(var, ae_dict, f_overlap)
            if ae is not None:
                if ae.startswith('SINE') or ae.startswith('LINE') or ae.split('|')[2].startswith('SVA'):
                    ae = 'ME:' + ae
                var.alt = '<DEL:%s>' % ae
                var.info['SVTYPE'] = 'MEI'
                vcf_out.write(var.get_var_string(True) + '\n')
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
                    #vcf_out.write(var.get_var_string(True) + '\n')
                    vcf_out.write(line)
                else:
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/low_freq_no_rd')
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/low_freq_no_rd')
                    for m_var in to_bnd_strings(var):
                        vcf_out.write(m_var + '\n')
            else:
                if has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold):
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/high_freq_rd')
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/high_freq_rd')
                    # write variant
                    #vcf_out.write(var.get_var_string(True) + '\n')
                    vcf_out.write(line)
                else:
                    # has_high_freq_depth_support(var, gender, exclude, slope_threshold, rsquared_threshold, writedir + '/high_freq_no_rd')
                    # has_low_freq_depth_support(var, gender, exclude, writedir + '/high_freq_no_rd')
                    for m_var in to_bnd_strings(var):
                        vcf_out.write(m_var + '\n')
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
    if args.exclude is not None:
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
