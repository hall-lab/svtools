#!/usr/bin/env python

import argparse, sys, copy, gzip, time, math, re
import numpy as np
import pandas as pd
from scipy import stats, cluster, spatial
from scipy.stats import multivariate_normal
from sklearn import metrics
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


class gauss_mix:
    
    def __init__(self, mns=None, covs=None, pp=None):
        self.__mu=mns
        self.__sigma=covs
        self.__probs=pp

    def __str__(self):
        string=str(self.__mu)+"\n"+str(self.__sigma)+"\n"+str(self.__probs)+"\n"
        return string

    def get_mu(self, i):
        return self.__mu[i]

    def set_mu(self, i, newmu_i):
        self.__mu[i]=newmu_i

    def get_sigma(self, i):
        return self.__sigma[i]

    def set_sigma(self, i, newsigma_i):
        self.__sigma[i]=newsigma_i

    def get_probs(self, i):
        return self.__probs[i]

    def set_probs(self, i, newpp_i):
        self.__probs[i]=newpp_i


    

def recluster(df):

    df=df[(df.AB!=".")].copy()
    df.loc[:,'AB']=pd.to_numeric(df.loc[:,'AB'])
    df.loc[:,'CN']=pd.to_numeric(df.loc[:,'CN'])

    tp=df.iloc[0,:].loc['svtype']
    gt_code={'0/0':1, '0/1':2, '1/1':3}
    gt_code_rev={1:'0/0', 2:'0/1', 3:'1/1'}
    df.loc[:,'gtn']=df.loc[:, 'GT'].map(gt_code)

    if tp in ['DEL']:
        recluster_DEL(df)
        df.loc[:,'GTR']=df.loc[:, 'gt_new'].map(gt_code_rev)
    elif tp in ['DUP']:
        recluster_DUP(df)
        re_recluster_DUP(df)
        df.loc[:,'GTR']=df.loc[:, 'GT'].copy()
    return df


def recluster_DEL(df):

    #priors
    mu_0={1:np.array([0.03, 2]),  2:np.array([0.32, 1.11]), 3:np.array([0.88, 0.20])}
    psi={ 1:np.matrix('0.02 -0.0001; -0.0001 0.12'),
      2:np.matrix('0.05 -0.02; -0.02 0.2'),
      3:np.matrix('0.15 -0.02; -0.02 0.1')}

    lambda_0=1
    nu_0=1

    gpd=df.loc[:, ['gtn', 'CN', 'AB']].groupby(['gtn'])
    covs=gpd[['AB','CN']].cov()
    mns=gpd[['AB', 'CN']].mean()
    cts=gpd.size()

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        lin_fit=smf.ols('CN~AB',df).fit()
        df.loc[:, 'gt_adj']=df.loc[:, 'gtn'].copy()
        #check that CN, AB are correlated, and in the right direction
        if (lin_fit.rsquared>0.5) and (-1*lin_fit.params[1]>0.5):
             x_int=-lin_fit.params[0]/lin_fit.params[1]
             #adjust init GT calls if AB shifted toward 0 
             if x_int<1:
                 #find mdpts between neighboring GT
                 mins=gpd['AB'].min()
                 maxes=gpd['AB'].max()
                 bound1=0.2 
                 bound2=0.7 
                 if (2 in mins) and (1 in maxes):
                     bound1=0.5*(mins[2]+maxes[1])
                 if (3 in mins) and (2 in maxes):
                     bound2=0.5*(mins[3]+maxes[2])
                 newbound1=bound1*x_int
                 newbound2=bound2*x_int
                 df.loc[:, 'gt_adj']=pd.to_numeric(pd.cut(df['AB'], bins=[-1, newbound1, newbound2, 1], labels=['1', '2', '3']))
                 gpd=df.loc[:,['gt_adj', 'CN', 'AB']].groupby(['gt_adj'])
                 covs=gpd[['AB', 'CN']].cov()
                 mns=gpd[['AB', 'CN']].mean()
                 cts=gpd.size()

    mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
            2: get_mu_map(2, cts, lambda_0, mu_0, mns),
            3: get_mu_map(3, cts, lambda_0, mu_0, mns)}

    sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
               2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
               3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}

    gm=gauss_mix(mu_map, sigma_map, {1:0.34, 2:0.33, 3:0.33})
    
    Estep(df, gm)
    Mstep(df, gm)
    old_ll=log_lik(df, gm)

    max_its=100
    tol=0.1

    for i in range(max_its):
        #sys.stderr.write(str(i)+"\t"+str(gm)+"\n")
        Estep(df, gm)
        Mstep(df, gm)
        ll=log_lik(df, gm)
        sys.stderr.write(str(i)+"\t"+str(ll-old_ll)+"\n")
        if(abs(ll-old_ll) < tol):
            sys.stderr.write("done\n")
            break
        old_ll=ll

    lld_code={'lld1':1, 'lld2':2, 'lld3':3}
    df.loc[:,'gt_new']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
    df.loc[:, 'gq']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
    df.loc[:, 'med_gq']=df.loc[:, 'gq'].median()    
    df.loc[:, 'q10_gq']=df.loc[:, 'gq'].quantile(0.1)

    return

def Estep(df1, gm1):

    df1.loc[:, 'lld1']=multivariate_normal.pdf(df1.loc[:, ['AB', 'CN']], mean=gm1.get_mu(1), cov=gm1.get_sigma(1), allow_singular=False)
    df1.loc[:, 'lld2']=multivariate_normal.pdf(df1.loc[:, ['AB', 'CN']], mean=gm1.get_mu(2), cov=gm1.get_sigma(2), allow_singular=False)
    df1.loc[:, 'lld3']=multivariate_normal.pdf(df1.loc[:, ['AB', 'CN']], mean=gm1.get_mu(3), cov=gm1.get_sigma(3), allow_singular=False)
    denom=df1.loc[:,'lld1']*gm1.get_probs(1)+df1.loc[:, 'lld2']*gm1.get_probs(2)+df1.loc[:, 'lld3']*gm1.get_probs(3)
    df1.loc[:, 'p1']=df1.loc[:, 'lld1']*gm1.get_probs(1)/denom
    df1.loc[:, 'p2']=df1.loc[:, 'lld2']*gm1.get_probs(2)/denom
    df1.loc[:, 'p3']=df1.loc[:, 'lld3']*gm1.get_probs(3)/denom
    
    return

def Mstep(df1, gm1):

    idmat=np.matrix('1 0; 0 1')
    eps=1e-3
    ridge=idmat*eps

    ss1=sws.DescrStatsW(df1[ ['AB', 'CN']], weights=df1.loc[:, 'p1'])
    ss2=sws.DescrStatsW(df1[['AB', 'CN']], weights=df1.loc[:, 'p2'])
    ss3=sws.DescrStatsW(df1[['AB', 'CN']], weights=df1.loc[:, 'p3'])

    gm1.set_mu(1, ss1.mean)
    gm1.set_mu(2, ss2.mean)
    gm1.set_mu(3, ss3.mean)
    gm1.set_sigma(1, np.add(ss1.cov, ridge))
    gm1.set_sigma(2, np.add(ss2.cov, ridge))
    gm1.set_sigma(3, np.add(ss3.cov, ridge))
    gm1.set_probs(1, df1['p1'].mean())
    gm1.set_probs(2, df1['p2'].mean())
    gm1.set_probs(3, df1['p3'].mean())


def log_lik(df1, gm1):
    df1.loc[:, 'lld']=(df1.loc[:, 'lld1']*gm1.get_probs(1)+df1.loc[:, 'lld2']*gm1.get_probs(2)+df1.loc[:, 'lld3']+gm1.get_probs(3))
    lld=df1.loc[:, 'lld'].sum()
    return lld


def re_recluster_DUP(df):

    df.loc[:, 'gt_new_re']=0
    df.loc[:, 'gq_re']=0
    df.loc[:, 'med_gq_re']=0
    df.loc[:, 'q10_gq_re']=0
    return

def re_recluster_DUP(df):

    #priors
    mu_0={1: np.array([0.03, 2]), 2:np.array([0.27,3]), 3:np.array([0.45,4])}
    psi={1:np.matrix('0.00128 -0.00075; -0.00075 1.1367'), 
      2:np.matrix('0.013 -0.0196; -0.0196 0.4626'),
      3:np.matrix('0.0046 -0.0112; -0.0112 0.07556')}
    lambda_0=1
    nu_0=1

    df.loc[:, 'gt_adj']=df.loc[:,'gt_new'].copy()
    df.loc[ (df['gt_new']==1) & (df['AB']>0.1) & (df['CN']>2.5), 'gt_adj']=2

    gpd=df.loc[:, ['gt_adj', 'CN', 'AB']].groupby(['gt_adj'])
    covs=gpd[['AB','CN']].cov()
    mns=gpd[['AB', 'CN']].mean()
    cts=gpd.size()

    mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
            2: get_mu_map(2, cts, lambda_0, mu_0, mns),
            3: get_mu_map(3, cts, lambda_0, mu_0, mns)}
    sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
               2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
               3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}

    df.loc[:, 'lld1']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[1], cov=sigma_map[1])
    df.loc[:, 'lld2']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[2], cov=sigma_map[2])
    df.loc[:, 'lld3']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[3], cov=sigma_map[3])
    lld_code={'lld1':1, 'lld2':2, 'lld3':3}
    df.loc[:,'gt_new_re']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
    df.loc[:, 'gq_re']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
    df.loc[:, 'med_gq_re']=df.loc[:, 'gq_re'].median()
    df.loc[:, 'q10_gq_re']=df.loc[:, 'gq_re'].quantile(0.1)
    return 

def recluster_DUP(df):

    #priors
    mu_0={1: np.array([0.03, 2]), 2:np.array([0.27,3]), 3:np.array([0.45,4])}
    psi={1:np.matrix('0.00128 -0.00075; -0.00075 1.1367'), 
      2:np.matrix('0.013 -0.0196; -0.0196 0.4626'),
      3:np.matrix('0.0046 -0.0112; -0.0112 0.07556')}
    lambda_0=1
    nu_0=1

    gpd=df.loc[:, ['gtn', 'CN', 'AB']].groupby(['gtn'])
    covs=gpd[['AB','CN']].cov()
    mns=gpd[['AB', 'CN']].mean()
    print(str(mns))
    exit(1)
    cts=gpd.size()
    
    df.loc[:, 'gt_adj']=df.loc[:, 'gtn'].copy()

    mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
            2: get_mu_map(2, cts, lambda_0, mu_0, mns),
            3: get_mu_map(3, cts, lambda_0, mu_0, mns)}
    sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
               2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
               3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}

    df.loc[:, 'lld1']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[1], cov=sigma_map[1])
    df.loc[:, 'lld2']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[2], cov=sigma_map[2])
    df.loc[:, 'lld3']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[3], cov=sigma_map[3])
    lld_code={'lld1':1, 'lld2':2, 'lld3':3}
    df.loc[:,'gt_new']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
    df.loc[:, 'gq']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
    df.loc[:, 'med_gq']=df.loc[:, 'gq'].median()
    df.loc[:, 'q10_gq']=df.loc[:, 'gq'].quantile(0.1)
    return 

def recluster_INV_BND(df):

    #priors
    mu_0={1: 0.03, 2:0.46, 3:0.94}
    psi={1:0.00128, 2:0.013, 3:0.0046}
    
    lambda_0=1
    nu_0=1

    gpd=df.loc[:, ['gtn', 'AB']].groupby(['gtn'])
    covs=gpd[['AB']].cov()
    mns=gpd[['AB']].mean()
    cts=gpd.size()
    
    df.loc[:, 'gt_adj']=df.loc[:, 'gtn'].copy()
    mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
            2: get_mu_map(2, cts, lambda_0, mu_0, mns),
            3: get_mu_map(3, cts, lambda_0, mu_0, mns)}
    sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
               2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
               3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}

    df.loc[:, 'lld1']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[1], cov=sigma_map[1])
    df.loc[:, 'lld2']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[2], cov=sigma_map[2])
    df.loc[:, 'lld3']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[3], cov=sigma_map[3])
    lld_code={'lld1':1, 'lld2':2, 'lld3':3}
    df.loc[:,'gt_new']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
    df.loc[:, 'gq']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
    df.loc[:, 'med_gq']=df.loc[:, 'gq'].median()
    df.loc[:, 'q10_gq']=df.loc[:, 'gq'].quantile(0.1)
    return 



def get_mu_map(j, nn, l0, m0, sample_mean):
    new_mu=m0[j]
    if j in nn:
        new_mu=(l0*m0[j]+nn[j]*sample_mean.loc[j,:].values)/(l0+nn[j])
    return new_mu

def get_sigma_map(j, nn, l0, psi0, sample_cov, sample_mean, m0):
    new_sig=psi0[j]
    if (j in nn) and nn[j]>2:
        val=sample_cov.loc[j,:]
        val=val+l0*nn[j]/(l0+nn[j])*np.outer(sample_mean.loc[j,:]-m0[j], sample_mean.loc[j,:]-m0[j])
        new_sig=val+psi0[j]
    return new_sig

        
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
                vcf.add_info('Q10GQR', '1', 'Float', 'Q10 quality for refined GT')
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
        # bail if not DEL or DUP prior to reclassification
        if svtype not in ['DEL']:
            vcf_out.write(line)
            continue
        
        var = Variant(v, vcf)
        sys.stderr.write("%s\n" % var.var_id)
        
        sys.stderr.write("%f\n" % float(var.get_info('AF')))
        if float(var.get_info('AF'))<0.01:
            vcf_out.write(line)
        else:
            df=load_df(var, exclude, sex)
            recdf=recluster(df)
            if ct==1:
                recdf.to_csv(outf, header=True)
                ct += 1
            else:
              recdf.to_csv(outf, header=False)
            #var.set_info("MEDGQR", '{:.2f}'.format(recdf.iloc[0,:].loc['med_gq_re']))
            #var.set_info("Q10GQR", '{:.2f}'.format(recdf.iloc[0,:].loc['q10_gq_re']))
            #recdf.set_index('sample', inplace=True)
            #for s in var.sample_list:
            #   if s in recdf.index:
            #        var.genotype(s).set_format("GTR", recdf.loc[s,'GTR'])
            #        var.genotype(s).set_format("GQR", '{:.2f}'.format(recdf.loc[s,'gq_re']))
            #    else:
            #        var.genotype(s).set_format("GTR", "./.")
            #        var.genotype(s).set_format("GQR", 0)
            #vcf_out.write(var.get_var_string(use_cached_gt_string=False) + '\n')

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
