#!/usr/bin/env python

import argparse, sys, copy, gzip, time, math, re
import numpy as np
import pandas as pd
from scipy import stats
from collections import Counter, defaultdict, namedtuple
import statsmodels.formula.api as smf
from operator import itemgetter
import warnings
from svtools.vcf.file import Vcf
from svtools.vcf.genotype import Genotype
from svtools.vcf.variant import Variant
import svtools.utils as su



def todosage(gt):
    if gt == '0/0':
        return 0
    elif gt == '0/1':
        return 1
    elif gt == '1/1':
        return 2
    else:
        return np.nan



def ld_calc(curlist, keep, ld_outfile, winsz):

    varlist=[var.var_id for var in curlist]
    df=pd.DataFrame(index=keep, columns=varlist)
    for var in curlist:
       gtlist=[todosage(var.genotype(s).get_format('GT')) for s in keep]
       df[var.var_id]=gtlist

    df.to_csv('./df.csv')
    cc=df.corr().stack().reset_index()
    cc.columns=['var1', 'var2', 'R']
    
    vars=list(df.columns.values)[1:]
    id=pd.DataFrame(vars, columns=['var'])
    id['ind']=id.index
    ccm=cc.merge(id, left_on='var1', right_on='var').merge(id, left_on='var2', right_on='var')[['var1', 'var2', 'R', 'ind_x', 'ind_y']]
    ccm['diff']=ccm['ind_y']-ccm['ind_x']
    ccm=ccm.loc[ccm['diff'].isin(range(winsz+1))].sort_values(by=['ind_x', 'ind_y'])
    ccm.to_csv('./ccm.csv', mode='a')




# primary function
def calc_ld(vcf_in,  exclude_file, ld_outfile, winsz, minpos):

    vcf = Vcf()
    header = []
    in_header = True
    maxwin=100

    exclude = []
    keep = []

    if exclude_file is not None:
        for line in exclude_file:
            exclude.append(line.rstrip())

    if ld_outfile is not None:
        outf=open(ld_outfile, 'w', 4096)
        outf.write("id1\tid2\tnp1\tnp2\tr2\n")

    curlist = []
    curchr = -1
    
    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                continue
            else:
                in_header = False
                vcf.add_header(header)

        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        for s in var.sample_list:
            if s in exclude:
                continue
            keep.append(s)

        if var.info['NSAMP']>minpos:
            if curchr!=-1 and var.chr is not curchr: 
                ld_calc(curlist, keep, ld_outfile, winsz)
                curlist=[var]
                curchr=var.chr
            elif len(curlist)>maxwin:
                ld_calc(curlist, keep, ld_outfile, winsz)
                curlist=curlist[(maxwin-1-winsz):]
                curlist.append(var)
            else:
                curlist.append(var)


    ld_calc(curlist, keep, ld_outfile, winsz)
    if ld_outfile is not None:
        outf.close()
    vcf_in.close()
    if exclude_file is not None:
        exclude_file.close()

    return


def run_pairwise_ld(vcf_file,  exclude_list,  ld_outfile, winsz):
    calc_ld(vcf_file, exclude_list, ld_outfile, winsz, 10)

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', default=None, help='VCF input, sorted by chr, pos')
    parser.add_argument('-e', '--exclude', metavar='<FILE>', dest='exclude', type=argparse.FileType('r'), required=False, default=None, help='list of samples to exclude from classification algorithms')
    parser.add_argument('-w', '--winsz', metavar='<INT>', dest='winsz', type=int, default=100, help='window size (in variants)')
    parser.add_argument('-l', '--ld_outfile', metavar='<STRING>', dest='ld_outfile', type=str, default=None, required=True, help='ld output file')
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'calculate pairwise ld'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        run_pairwise_ld(stream, args.exclude,  args.ld_outfile, args.winsz)


if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
