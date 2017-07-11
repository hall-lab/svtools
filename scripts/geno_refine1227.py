#!/usr/bin/env python

import argparse, sys, copy, gzip, time, math, re
import numpy as np
import pandas as pd
from scipy import stats, cluster, spatial
from scipy.stats import multivariate_normal, norm
from sklearn import metrics, mixture
from collections import Counter, defaultdict, namedtuple
import statsmodels.formula.api as smf
import statsmodels.stats.weightstats as sws
from operator import itemgetter
import warnings
from svtools.vcf.file import Vcf
from svtools.vcf.genotype import Genotype
from svtools.vcf.variant import Variant
import svtools.utils as su

vcf_rec = namedtuple ('vcf_rec', 'var_id sample svtype AF GT CN AB ')
    

def recluster(df):

    df=df[(df.AB!=".")].copy()
    df.loc[:,'AB']=pd.to_numeric(df.loc[:,'AB'])
    df.loc[:,'CN']=pd.to_numeric(df.loc[:,'CN'])

    tp=df.iloc[0,:].loc['svtype']
    gt_code={'./.':0, '0/0':1, '0/1':2, '1/1':3}
    df.loc[:,'gtn']=df.loc[:, 'GT'].map(gt_code)

    if tp in ['DEL', 'MEI']:
        recluster_DEL(df)
    elif tp in ['DUP']:
        recluster_DUP(df)
    elif tp in ['INV']:
        recluster_INV_BND(df)
    elif tp in ['BND']:
        recluster_INV_BND(df)
    return df


def recluster_DEL(df):

    mu0=np.array([[0.03, 2.],
                 [0.32, 1.11],
                 [0.88, 0.20]])
    sigma_inv0=np.array([[50.0, 5.0], 
                         [20.0, 5.0],
                         [6.7, 5.0]])
    dat=df.loc[:, [ 'AB', 'CN']].as_matrix()
    cts=np.bincount(df['gtn'], minlength=4)                    
    pp=math.sqrt(1.0*cts[1]/cts.sum())
    
    if pp<0.05:  #just punt for now
        df['gq']=0
        df['gt_new']=df['gtn'].copy()
        df['med_gq']=0
        df['q20_gq']=0
        df['GTR']='./.'
        return 

    wts0=np.array([pp*pp, 2*pp*(1-pp), (1-pp)*(1-pp)])+0.01
    gmm2=mixture.GaussianMixture(n_components=2, covariance_type='diag', reg_covar=0.0001, weights_init=wts0[:2]/wts0[:2].sum(), means_init=mu0[:2], precisions_init=sigma_inv0[:2], random_state=1 ).fit(dat)
    gmm3=mixture.GaussianMixture(n_components=3, covariance_type='diag', reg_covar=0.0001, weights_init=wts0/wts0.sum(), means_init=mu0, precisions_init=sigma_inv0, random_state=1 ).fit(dat)
    gmm1=mixture.GaussianMixture(n_components=1, covariance_type='diag', reg_covar=0.0001, random_state=1).fit(dat)
    bic=np.array([gmm1.bic(dat), gmm2.bic(dat), gmm3.bic(dat)])
    gmms=np.array([gmm1, gmm2, gmm3])

    model=gmms[bic.argmin()]
    gq=calc_gq(dat, model)
    gt=model.predict(dat)
    df['gq']=gq
    df['gt_new']=gt
    df['med_gq']=np.percentile(gq, 50)
    df['q20_gq']=np.percentile(gq, 20)                                                                                                                                                         
    cn_means=model.means_[:,1]
    ordr=np.argsort(cn_means)[::-1]
    if cn_means.size==1:
        ordr=np.append(ordr, 1)
        cn_means=np.append(cn_means, 1)
    if cn_means.size==2:
        ordr=np.append(ordr, 2)
        cn_means=np.append(cn_means, 0)

    gt_code_rev={ordr[0]:'0/0', ordr[1]:'0/1', ordr[2]:'1/1'}
    df.loc[:,'GTR']=df.loc[:, 'gt_new'].map(gt_code_rev)
  
    return    

def recluster_DUP(df):

    dat=df.loc[:, [ 'AB', 'CN']].as_matrix()
    dat1=dat[:,1].reshape(-1, 1)
    model_multi=recluster_DUP_multi(df)
    gq_multi=calc_gq(dat1, model_multi)
    gq_multi_med=np.percentile(gq_multi, 50)

    model_bi=recluster_DUP_bi(df)
    gq_bi=calc_gq(dat, model_bi)
    gq_bi_med=np.percentile(gq_bi, 50)

    if gq_bi_med>gq_multi_med:
        gt=model_bi.predict(dat)
        gq=gq_bi
        gq_med=gq_bi_med
    else:
        gt=model_multi.predict(dat1)
        gq=gq_multi
        gq_med=gq_multi_med

    df['gq']=gq
    df['gt_new']=gt
    df['med_gq']=gq_med
    df['q20_gq']=np.percentile(gq, 20)

    genotype=False
    if np.percentile(dat[:,1], 99)<4.5:
        
        if gq_bi_med>gq_multi_med:
            cn_means=model_bi.means_[:,1]
        else:
            cn_means=model_multi.means_[:,0]
        if cn_means.size<=3:
            ordr=np.argsort(cn_means)
            if cn_means.size==1:
                ordr=np.append(ordr, 1)
                cn_means=np.append(cn_means,3)
            if cn_means.size==2:
                ordr=np.append(ordr, 2)
                cn_means=np.append(cn_means, 4)
            if abs(cn_means[ordr[0]]-2)<1 and abs(cn_means[ordr[1]]-3)<1 and  abs(cn_means[ordr[2]]-4)<1.5:
                genotype=True
                gt_code_rev={ordr[0]:'0/0', ordr[1]:'0/1', ordr[2]:'1/1'}
                df.loc[:,'GTR']=df.loc[:, 'gt_new'].map(gt_code_rev)
    if genotype is False:
        df.loc[:, 'GTR']='./.'

    return
            

def recluster_DUP_multi(df):
    
    dat=df.loc[:, [ 'AB', 'CN']].as_matrix()
    dat1=dat[:,1].reshape(-1, 1)
    gmm=[]
    bic=[]
    for i in range(1,6):
        model=mixture.GaussianMixture(n_components=i, covariance_type='diag', reg_covar=0.0001).fit(dat1)
        gmm.append(model)
        bic.append(model.bic(dat1))

    bic=np.array(bic)
    model=gmm[bic.argmin()]
    return model

def recluster_DUP_bi(df):

    mu0=np.array([[0.03, 2.],
                 [0.2, 2.7],
                 [0.6, 3.5]])

    sigma_inv0=np.array([[50.0, 5.0], 
                         [20.0, 5.0],
                         [6.7, 5.0]])

    dat=df.loc[:, [ 'AB', 'CN']].as_matrix()
    cts=np.bincount(df['gtn'], minlength=4)                    
    
    pp=math.sqrt(1.0*cts[1]/cts.sum())
  
    wts0=np.array([pp*pp, 2*pp*(1-pp), (1-pp)*(1-pp)])+0.01
    gmm2=mixture.GaussianMixture(n_components=2, covariance_type='diag', reg_covar=0.0001, weights_init=wts0[:2]/wts0[:2].sum(), means_init=mu0[:2], precisions_init=sigma_inv0[:2], random_state=1 ).fit(dat)
    gmm3=mixture.GaussianMixture(n_components=3, covariance_type='diag', reg_covar=0.0001, weights_init=wts0/wts0.sum(), means_init=mu0, precisions_init=sigma_inv0, random_state=1 ).fit(dat)
    gmm1=mixture.GaussianMixture(n_components=1, covariance_type='diag', reg_covar=0.0001, random_state=1).fit(dat)
    bic=np.array([gmm1.bic(dat), gmm2.bic(dat), gmm3.bic(dat)])
    gmms=np.array([gmm1, gmm2, gmm3])

    model=gmms[bic.argmin()]
    return model


def calc_gq(dat, model):

    probs=model.predict_proba(dat)
    nclus=probs.shape[1]
    ranks=np.argsort(probs, axis=1)
    wbest=model.weights_[ranks[:,(nclus-1)]]
    w2nd=model.weights_[ranks[:,(nclus-2)]]
    pbest=probs[np.arange(len(probs)), ranks[:,(nclus-1)]]
    p2nd=probs[np.arange(len(probs)), ranks[:,(nclus-2)]]
    epsilon=1e-200
    gq=np.log10(pbest+epsilon)-np.log10(wbest+epsilon)-(np.log10(p2nd+epsilon)-np.log10(w2nd+epsilon))
    return gq

def recluster_INV_BND(df):

    dim=1
    mu0=np.array([[0.3], [0.3], [0.7]])
    sigma_inv0=np.array([[50.0], [20.0], [6.7]])
    dat=df.loc[:, [ 'AB']].as_matrix().reshape(-1, 1)
    
    cts=np.bincount(df['gtn'], minlength=4)

    pp=math.sqrt(1.0*cts[1]/cts.sum())
    if pp<0.05:  #just punt for now                                                                                                                                                                                                              
        df['gq']=0
        df['gt_new']=df['gtn'].copy()
        df['med_gq']=0
        df['q20_gq']=0
        df['GTR']='./.'
        return

    wts0=np.array([pp*pp, 2*pp*(1-pp), (1-pp)*(1-pp)])+0.01
    gmm2=mixture.GaussianMixture(n_components=2, covariance_type='diag', reg_covar=0.0001, weights_init=wts0[:2]/wts0[:2].sum(), means_init=mu0[:2], precisions_init=sigma_inv0[:2], random_state=1 ).fit(dat)
    gmm3=mixture.GaussianMixture(n_components=3, covariance_type='diag', reg_covar=0.0001, weights_init=wts0/wts0.sum(), means_init=mu0, precisions_init=sigma_inv0, random_state=1 ).fit(dat)
    gmm1=mixture.GaussianMixture(n_components=1, covariance_type='diag', reg_covar=0.0001, random_state=1).fit(dat)
    bic=np.array([gmm1.bic(dat), gmm2.bic(dat), gmm3.bic(dat)])
    gmms=np.array([gmm1, gmm2, gmm3])

    model=gmms[bic.argmin()]
    gq=calc_gq(dat, model)
    gt=model.predict(dat)
    df['gq']=gq
    df['gt_new']=gt
    df['med_gq']=np.percentile(gq, 50)
    df['q20_gq']=np.percentile(gq, 20)
    ab_means=model.means_[:,0]
    ordr=np.argsort(ab_means)
    
    if ab_means.size==1:
        ordr=np.append(ordr, 1)
        ab_means=np.append(ab_means, 1)
    if ab_means.size==2:
        ordr=np.append(ordr, 2)
        ab_mean=np.append(ab_means, 1)
    gt_code_rev={ordr[0]:'0/0', ordr[1]:'0/1', ordr[2]:'1/1'}
    df.loc[:,'GTR']=df.loc[:, 'gt_new'].map(gt_code_rev)
       
    return

#def recluster_INV_BND(df):
#
#
#    if(np.cov(df['CN'])>0):
#        dim=2
#        mu0=np.array([[0.03, 2.],
#                 [0.3, 2.],
#                 [0.88, 2.]])
#
#        sigma_inv0=np.array([[50.0, 5.0], 
#                         [20.0, 5.0],
#                         [6.7, 5.0]])
#        dat=df.loc[:, [ 'AB', 'CN']].as_matrix()
#
#    else:
#        dim=1
#        mu0=np.array([[0.3], [0.3], [0.7]])
#        sigma_inv0=np.array([[50.0], [20.0], [6.7]])
#        dat=df.loc[:, [ 'AB']].as_matrix().reshape(-1, 1)
#    
#    cts=np.bincount(df['gtn'], minlength=4)
#
#    pp=math.sqrt(1.0*cts[1]/cts.sum())
#    if pp<0.05:  #just punt for now                                                                                                                                                                                                              
#        df['gq']=0
#        df['gt_new']=df['gtn'].copy()
#        df['med_gq']=0
#        df['q20_gq']=0
#        df['GTR']='./.'
#        return
#
#    wts0=np.array([pp*pp, 2*pp*(1-pp), (1-pp)*(1-pp)])+0.01
#    gmm2=mixture.GaussianMixture(n_components=2, covariance_type='diag', reg_covar=0.0001, weights_init=wts0[:2]/wts0[:2].sum(), means_init=mu0[:2], precisions_init=sigma_inv0[:2], random_state=1 ).fit(dat)
#    gmm3=mixture.GaussianMixture(n_components=3, covariance_type='diag', reg_covar=0.0001, weights_init=wts0/wts0.sum(), means_init=mu0, precisions_init=sigma_inv0, random_state=1 ).fit(dat)
#    gmm1=mixture.GaussianMixture(n_components=1, covariance_type='diag', reg_covar=0.0001, random_state=1).fit(dat)
#    bic=np.array([gmm1.bic(dat), gmm2.bic(dat), gmm3.bic(dat)])
#    gmms=np.array([gmm1, gmm2, gmm3])
#
#    model=gmms[bic.argmin()]
#    gq=calc_gq(dat, model)
#    gt=model.predict(dat)
#    df['gq']=gq
#    df['gt_new']=gt
#    df['med_gq']=np.percentile(gq, 50)
#    df['q20_gq']=np.percentile(gq, 20)
#    ab_means=model.means_[:,0]
#    ordr=np.argsort(ab_means)
#    
#    if ab_means[ordr[0]]<0.01:
#        if ab_means.size==1:
#            ordr=np.append(ordr, 1)
#            ab_means=np.append(ab_means, 1)
#        if ab_means.size==2:
#            ordr=np.append(ordr, 2)
#            ab_mean=np.append(ab_means, 1)
#        gt_code_rev={ordr[0]:'0/0', ordr[1]:'0/1', ordr[2]:'1/1'}
#        df.loc[:,'GTR']=df.loc[:, 'gt_new'].map(gt_code_rev)
#    else:
#        df.loc[:, 'GTR']='./.'
#       
#    return


        
def percentile(n):
    def percentile_(x):
        return np.percentile(x,n)
    percentile_.__name__= 'percentile_%s' % n
    return percentile_


def load_df(var, exclude, sex):
    test_set = list()
    for s in var.sample_list:
        if s in exclude:
            continue
        if 'CN' not in var.format_dict:
            cn=0
        else:
            cn = var.genotype(s).get_format('CN')
            if (var.chrom == 'X' or var.chrom == 'Y') and sex[s] == 1:
                cn=str(float(cn)*2)
        test_set.append(vcf_rec(var.var_id, s, var.info['SVTYPE'], var.info['AF'], 
                    var.genotype(s).get_format('GT'),  cn , var.genotype(s).get_format('AB')))

    test_set = pd.DataFrame(data = test_set, columns=vcf_rec._fields)
    return test_set


def run_gt_refine(vcf_in, vcf_out, diag_outfile, gender_file, exclude_file):

    vcf = Vcf()
    header = []
    in_header = True
    sex={}
    
    for line in gender_file:
        v = line.rstrip().split('\t')
        sex[v[0]] = int(v[1])

    exclude = []
    if exclude_file is not None:
        for line in exclude_file:
            exclude.append(line.rstrip())

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
                vcf.add_info('MEDGQR', '1', 'Float', 'Median quality for refined GT')
                vcf.add_info('Q20GQR', '1', 'Float', 'Q20 quality for refined GT')
                vcf.add_format('GQR', 1, 'Float', 'Quality of refined genotype.')
                vcf.add_format('GTR', 1, 'String', 'Refined genotype.')
                vcf_out.write(vcf.get_header() + '\n')

        v = line.rstrip().split('\t')
        info = v[7].split(';')
        svtype = None
        for x in info:
            if x.startswith('SVTYPE='):
                svtype = x.split('=')[1]
                break
        
        var = Variant(v, vcf)
        sys.stderr.write("%s\n" % var.var_id)
        
        sys.stderr.write("%f\n" % float(var.get_info('AF')))
        df=load_df(var, exclude, sex)
        recdf=recluster(df)
        if ct==1:
            recdf.to_csv(outf, header=True)
            ct += 1
        else:
            recdf.to_csv(outf, header=False)
        var.set_info("MEDGQR", '{:.2f}'.format(recdf.iloc[0,:].loc['med_gq']))
        var.set_info("Q20GQR", '{:.2f}'.format(recdf.iloc[0,:].loc['q20_gq']))
        recdf.set_index('sample', inplace=True)
        for s in var.sample_list:
            if s in recdf.index:
                var.genotype(s).set_format("GTR", recdf.loc[s,'GTR'])
                var.genotype(s).set_format("GQR", '{:.2f}'.format(recdf.loc[s,'gq']))
            else:
                var.genotype(s).set_format("GTR", "./.")
                var.genotype(s).set_format("GQR", 0)
                
        vcf_out.write(var.get_var_string(use_cached_gt_string=False) + '\n')

    vcf_out.close()
    vcf_in.close()
    gender_file.close()
    outf.close()
    
    if exclude_file is not None:
        exclude_file.close()
    return

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
    parser.add_argument('-d', '--diag_file', metavar='<STRING>', dest='diag_outfile', type=str, default=None, required=False, help='text file to output method comparisons')
    parser.add_argument('-g', '--gender', metavar='<FILE>', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
    parser.add_argument('-e', '--exclude', metavar='<FILE>', dest='exclude', type=argparse.FileType('r'), required=False, default=None, help='list of samples to exclude from classification algorithms')
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'refine genotypes by clustering'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.vcf_in) as stream:
        run_gt_refine(stream, args.vcf_out, args.diag_outfile, args.gender, args.exclude)

if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
