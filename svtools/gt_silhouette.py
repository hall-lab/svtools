#!/usr/bin/env python

import argparse, sys, copy, gzip, time, math, re
import numpy as np
import pandas as pd
from scipy import stats, cluster, spatial
from sklearn import metrics
from collections import Counter, defaultdict, namedtuple
import statsmodels.formula.api as smf
from operator import itemgetter
import warnings
from svtools.vcf.file import Vcf
from svtools.vcf.genotype import Genotype
from svtools.vcf.variant import Variant
import svtools.utils as su

vcf_rec = namedtuple ('vcf_rec', 'var_id sample svtype AF GT CN AB')

def get_silhouette(df):
    df=df[(df.AB!=".")].copy()
    df.loc[:,'AB']=pd.to_numeric(df.loc[:,'AB'])
    df.loc[:,'CN']=pd.to_numeric(df.loc[:,'CN'])

    tp=df.iloc[0,:].loc['svtype']

    [mn_CN, mn_AB]=df.loc[:, ['CN', 'AB']].mean(skipna=True)
    [sd_CN, sd_AB]=df.loc[:, ['CN', 'AB']].std(skipna=True)

    #standardize the 2 dims
    if sd_AB>0.01:
        df.loc[:, 'AB1']=(df.loc[:,'AB']-mn_AB)/sd_AB
    else: 
        df.loc[:, 'AB1']=df.loc[:, 'AB']
    if tp in ['DEL', 'DUP', 'MEI'] or sd_CN>0.01:
        df.loc[:, 'CN1']=(df.loc[:,'CN']-mn_CN)/sd_CN
    else:
        df.loc[:, 'CN1']=df.loc[:, 'CN']

    
    gt_code={'0/0':1, '0/1':2, '1/1':3}
    df.loc[:,'gtn']=df.loc[:, 'GT'].map(gt_code)

    dist_2d_sq=spatial.distance.squareform(spatial.distance.pdist(df[['AB1', 'CN1']], metric='cityblock'))
    df.loc[:, 'sil_gt_avg']=metrics.silhouette_score(dist_2d_sq, df.loc[:, 'gtn'].values, metric='precomputed')
    df.loc[:, 'sil_gt']=metrics.silhouette_samples(dist_2d_sq, df.loc[:, 'gtn'].values, metric='precomputed')
    df=df[ ['var_id', 'sample', 'svtype', 'AF', 'GT', 'CN', 'AB', 'sil_gt_avg', 'sil_gt']]
    return df

def load_df(var, sex):
    test_set = list()
    for s in var.sample_list:
        svtype=var.get_info('SVTYPE')
        if svtype in ['DEL', 'DUP', 'MEI']:
            cn = var.genotype(s).get_format('CN')
        else:
            cn=0
        if (var.chrom == 'X' or var.chrom == 'Y') and sex[s] == 1:
            cn=str(float(cn)*2)
        test_set.append(vcf_rec(var.var_id, s, var.info['SVTYPE'], var.info['AF'],
             var.genotype(s).get_format('GT'),  cn , var.genotype(s).get_format('AB')))

    test_set = pd.DataFrame(data = test_set, columns=vcf_rec._fields)
    return test_set


def run_gt_refine(vcf_in, vcf_out, diag_outfile, gender_file):

    vcf = Vcf()
    header = []
    in_header = True
    sex={}

    for line in gender_file:
        v = line.rstrip().split('\t')
        sex[v[0]] = int(v[1])

    outf=open(diag_outfile, 'w', 4096)
    ct=1
    
    for line in vcf_in:
        if in_header:
            if line[0] == "#":
               header.append(line)
               continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf.add_info('SIL_GT_AVG', '1', 'Float', 'Average silhouette of genotype clusters')
                #vcf.add_format('SIL_GT', '1', 'Float', 'Per-sample genotype cluster silhouette')
                vcf_out.write(vcf.get_header() + '\n')

        var = Variant(line.rstrip().split('\t'), vcf)
        df=load_df(var,  sex)
        df1=get_silhouette(df)

        sil_avg=df1.iloc[0, df1.columns.get_loc('sil_gt_avg')]
        #sil_ind=df1.loc[:, 'sil_gt']
        var.info['SIL_GT_AVG'] = '%0.2f' % sil_avg
        vcf_out.write(var.get_var_string(use_cached_gt_string=True) + '\n')
        
        if ct==1:
            df1.to_csv(outf, header=True)
            ct += 1
        else:
            df1.to_csv(outf, header=False)

    vcf_out.close()
    vcf_in.close()
    outf.close()
    gender_file.close()

    return

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
    parser.add_argument('-d', '--diag_file', metavar='<STRING>', dest='diag_outfile', type=str, default=None, required=False, help='text file to output method comparisons')
    parser.add_argument('-g', '--gender', metavar='<FILE>', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'refine genotypes by clustering'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.vcf_in) as stream:
        run_gt_refine(stream, args.vcf_out, args.diag_outfile, args.gender)


if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
