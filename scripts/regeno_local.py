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
#from svtools.vcf.file import Vcf
#from svtools.vcf.genotype import Genotype
#from svtools.vcf.variant import Variant
#import svtools.utils as su

vcf_rec = namedtuple ('vcf_rec', 'var_id sample svtype AF GT CN AB ')


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
    elif tp in ['BND']:
        recluster_BND(df)
        df.loc[:, 'GTR']=df.loc[:, 'gt_new'].map(gt_code_rev)
    elif tp in ['DUP']:
        recluster_DUP(df)
        df.loc[:,'GTR']=df.loc[:, 'gt_new'].copy()
    return df


def recluster_DEL(df):


    mu0=np.array([[0.03, 2],
                 [0.32, 1.11],
                 [0.88, 0.20]])
    #sigma0=np.array([[0.02, 0.2],
    #                 [0.05, 0.2],
    #                 [0.15, 0.2]]

    sigma_inv0=np.array([[50.0, 5.0], 
                         [20.0, 5.0],
                         [6.7, 5.0]])
                         
    gpd=df.loc[:, ['gtn', 'CN', 'AB']].groupby(['gtn'])
    cts=gpd.size()
    pp=math.sqrt(cts[1]/cts.sum())
    wts0=np.array([pp*pp, 2*pp*(1-pp), (1-pp)*(1-pp)])


    dat=df.loc[:, ['CN', 'AB']].as_matrix()
    gmm2=mixture.GaussianMixture(n_components=2, covariance_type='diag', reg_covar=0.0001, 
                                 weights_init=wts0[:2]/wts0[:2].sum(), means_init=mu0[:2], precisions_init=sigma_inv0[:2], random_state=1 ).fit(dat)

    gmm3=mixture.GaussianMixture(n_components=3, covariance_type='diag', reg_covar=0.0001, 
                                 weights_init=wts0/wts0.sum(), means_init=mu0, precisions_init=sigma_inv0, random_state=1 ).fit(dat)

    gmm1=mixture.GaussianMixture(n_components=1, covariance_type='diag', reg_covar=0.0001,.fit(dat)
    bic=np.array([gmm1.bic(dat), gmm2.bic(dat), gmm3.bic(dat)])
    gmms=np.array([gmm1, gmm2, gmm3])
    model=gmms[bic.argmin()]

    probs=model.predict_proba(dat)
    nclus=probs.shape[1]
    ranks=np.argsort(probs, axis=1)
    wbest=wts[ranks[:,(nclus-1)]]
    w2nd=wts[ranks[:,(nclus-2)]]
    pbest=probs[np.arange(len(probs)), ranks[:,(nclus-1)]]
    p2nd=probs[np.arange(len(probs)), ranks[:,(nclus-2)]]
    epsilon=1e-200
    gq=np.log10(pbest+epsilon)-np.log10(wbest+epsilon)-(np.log10(p2nd+epsilon)-np.log10(w2nd+epsilon))
    gt=model.predict(dat)
    df['gq']=gq
    df['gt_new']=gt
    df['med_gq']=np.percentile(gq, 50)
    df['q10_gq']=np.percentile(gq, 20)                                                                                                                                                         

     
    return    


#def recluster_BND(df):
#    
#    eps=1e-4
#    df.loc[:, 'gt_adj']=df.loc[:, 'gtn'].copy()
#    gpd=df.loc[:,['gt_adj', 'AB']].groupby(['gt_adj'])
#    covs=gpd[['AB']].std()
#    mns=gpd[['AB']].mean()
#    cts=gpd.size()
#
#    ngps=len(cts)
#
#    if ngps==3:
#        mu_map={1: mns.loc[1,:],
#                2: mns.loc[2,:],
#                3: mns.loc[3,:]}
#        sigma_map={1: covs.loc[1,:]+eps,
#                   2: covs.loc[2,:]+eps,
#                   3: covs.loc[3,:]+eps}
#        gm=gauss_mix(mu_map, sigma_map, {1:0.34, 2:0.33, 3:0.33}, 3, 1)
#
#    elif ngps==2:
#        mu_map={1: mns.loc[1,:],
#                2: mns.loc[2:,]}
#        sigma_map={1: covs.loc[1,:]+eps,
#                   2: covs.loc[2,:]+eps}
#        gm=gauss_mix(mu_map, sigma_map, {1:0.5, 2:0.5}, 2, 1)
#    
#    sys.stderr.write(str(gm))
#    Estep(df, gm)
#    Mstep(df, gm)
#    old_ll=log_lik(df, gm)
#
#    max_its=100
#    tol=0.1
#
#    for i in range(max_its):
#
#        Estep(df, gm)
#        Mstep(df, gm)
#        ll=log_lik(df, gm)
#        sys.stderr.write(str(i)+"\t"+str(ll-old_ll)+"\n")
#        sys.stderr.write(str(gm)+"\n")
#        if(abs(ll-old_ll) < tol):
#            sys.stderr.write("done\n")
#            break
#        old_ll=ll
#
#    lld_code={'lld1':1, 'lld2':2, 'lld3':3}
#    df.loc[:,'gt_new']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
#    df.loc[:, 'gq']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
#    df.loc[:, 'med_gq']=df.loc[:, 'gq'].median()    
#    df.loc[:, 'q10_gq']=df.loc[:, 'gq'].quantile(0.1)
#
#    return
#
#
#ef calc_likelihoods(df1, gm1):
#   ncomp=gm1.get_ncomp()
#   dim=gm1.get_dim()
#   sys.stderr.write(str(dim)+"\n in calc\n")
#   if dim==2:
#       df1.loc[:, 'lld1']=multivariate_normal.pdf(df1.loc[:, ['AB', 'CN']], mean=gm1.get_mu(1), cov=gm1.get_sigma(1), allow_singular=False)
#       df1.loc[:, 'lld2']=multivariate_normal.pdf(df1.loc[:, ['AB', 'CN']], mean=gm1.get_mu(2), cov=gm1.get_sigma(2), allow_singular=False)
#       if ncomp==3:
#           df1.loc[:, 'lld3']=multivariate_normal.pdf(df1.loc[:, ['AB', 'CN']], mean=gm1.get_mu(3), cov=gm1.get_sigma(3), allow_singular=False)
#       else:
#           df1.loc[:, 'lld3']=0
#   elif dim==1:
#       sys.stderr.write("ok\n")
#       df1.loc[:, 'lld1']=norm.pdf(df1.loc[:, ['AB']], loc=gm1.get_mu(1), scale=gm1.get_sigma(1))
#       df1.loc[:, 'lld2']=norm.pdf(df1.loc[:, ['AB']], loc=gm1.get_mu(2), scale=gm1.get_sigma(2))
#       if ncomp==3:
#           df1.loc[:, 'lld3']=norm.pdf(df1.loc[:, ['AB']], loc=gm1.get_mu(3), scale=gm1.get_sigma(3))
#       else:
#           df1.loc[:, 'lld3']=0



#ef re_recluster_DUP(df):
#
#   df.loc[:, 'gt_new_re']=0
#   df.loc[:, 'gq_re']=0
#   df.loc[:, 'med_gq_re']=0
#   df.loc[:, 'q10_gq_re']=0
#   return
#
#
#ef re_recluster_DUP(df):
#
#   #priors
#   mu_0={1: np.array([0.03, 2]), 2:np.array([0.27,3]), 3:np.array([0.45,4])}
#   psi={1:np.matrix('0.00128 -0.00075; -0.00075 1.1367'), 
#     2:np.matrix('0.013 -0.0196; -0.0196 0.4626'),
#     3:np.matrix('0.0046 -0.0112; -0.0112 0.07556')}
#   lambda_0=1
#   nu_0=1
#
#   df.loc[:, 'gt_adj']=df.loc[:,'gt_new'].copy()
#   df.loc[ (df['gt_new']==1) & (df['AB']>0.1) & (df['CN']>2.5), 'gt_adj']=2
#
#   gpd=df.loc[:, ['gt_adj', 'CN', 'AB']].groupby(['gt_adj'])
#   covs=gpd[['AB','CN']].cov()
#   mns=gpd[['AB', 'CN']].mean()
#   cts=gpd.size()
#
#   mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
#           2: get_mu_map(2, cts, lambda_0, mu_0, mns),
#           3: get_mu_map(3, cts, lambda_0, mu_0, mns)}
#   sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
#              2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
#              3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}
#
#   df.loc[:, 'lld1']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[1], cov=sigma_map[1])
#   df.loc[:, 'lld2']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[2], cov=sigma_map[2])
#   df.loc[:, 'lld3']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[3], cov=sigma_map[3])
#   lld_code={'lld1':1, 'lld2':2, 'lld3':3}
#   df.loc[:,'gt_new_re']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
#   df.loc[:, 'gq_re']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
#   df.loc[:, 'med_gq_re']=df.loc[:, 'gq_re'].median()
#   df.loc[:, 'q10_gq_re']=df.loc[:, 'gq_re'].quantile(0.1)
#   return 
#
#
#ef recluster_DUP(df):
#
#   #priors
#   mu_0={1: np.array([0.03, 2]), 2:np.array([0.27,3]), 3:np.array([0.45,4])}
#   psi={1:np.matrix('0.00128 -0.00075; -0.00075 1.1367'), 
#     2:np.matrix('0.013 -0.0196; -0.0196 0.4626'),
#     3:np.matrix('0.0046 -0.0112; -0.0112 0.07556')}
#   lambda_0=1
#   nu_0=1
#
#   gpd=df.loc[:, ['gtn', 'CN', 'AB']].groupby(['gtn'])
#   covs=gpd[['AB','CN']].cov()
#   mns=gpd[['AB', 'CN']].mean()
#   print(str(mns))
#   exit(1)
#   cts=gpd.size()
#   
#   df.loc[:, 'gt_adj']=df.loc[:, 'gtn'].copy()
#
#   mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
#           2: get_mu_map(2, cts, lambda_0, mu_0, mns),
#           3: get_mu_map(3, cts, lambda_0, mu_0, mns)}
#   sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
#              2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
#              3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}
#
#   df.loc[:, 'lld1']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[1], cov=sigma_map[1])
#   df.loc[:, 'lld2']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[2], cov=sigma_map[2])
#   df.loc[:, 'lld3']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[3], cov=sigma_map[3])
#   lld_code={'lld1':1, 'lld2':2, 'lld3':3}
#   df.loc[:,'gt_new']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
#   df.loc[:, 'gq']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
#   df.loc[:, 'med_gq']=df.loc[:, 'gq'].median()
#   df.loc[:, 'q10_gq']=df.loc[:, 'gq'].quantile(0.1)
#   return 
#
#
#ef recluster_INV_BND(df):
#
#   #priors
#   mu_0={1: 0.03, 2:0.46, 3:0.94}
#   psi={1:0.00128, 2:0.013, 3:0.0046}
#   
#   lambda_0=1
#   nu_0=1
#
#   gpd=df.loc[:, ['gtn', 'AB']].groupby(['gtn'])
#   covs=gpd[['AB']].cov()
#   mns=gpd[['AB']].mean()
#   cts=gpd.size()
#   
#   df.loc[:, 'gt_adj']=df.loc[:, 'gtn'].copy()
#   mu_map={1: get_mu_map(1, cts, lambda_0, mu_0, mns),
#           2: get_mu_map(2, cts, lambda_0, mu_0, mns),
#           3: get_mu_map(3, cts, lambda_0, mu_0, mns)}
#   sigma_map={1: get_sigma_map(1, cts, lambda_0, psi, covs, mns, mu_0),
#              2: get_sigma_map(2, cts, lambda_0, psi, covs, mns, mu_0),
#              3: get_sigma_map(3, cts, lambda_0, psi, covs, mns, mu_0)}
#
#   df.loc[:, 'lld1']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[1], cov=sigma_map[1])
#   df.loc[:, 'lld2']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[2], cov=sigma_map[2])
#   df.loc[:, 'lld3']=multivariate_normal.logpdf(df.loc[:, ['AB', 'CN']], mean=mu_map[3], cov=sigma_map[3])
#   lld_code={'lld1':1, 'lld2':2, 'lld3':3}
#   df.loc[:,'gt_new']=df.loc[:, ['lld1', 'lld2', 'lld3']].idxmax(1).map(lld_code)
#   df.loc[:, 'gq']=df.loc[:, ['lld1', 'lld2', 'lld3']].max(axis=1)-df.loc[:, ['lld1', 'lld2', 'lld3']].median(axis=1)
#   df.loc[:, 'med_gq']=df.loc[:, 'gq'].median()
#   df.loc[:, 'q10_gq']=df.loc[:, 'gq'].quantile(0.1)
#   return 
#
#
#ef get_mu_map(j, nn, l0, m0, sample_mean):
#   new_mu=m0[j]
#   if j in nn:
#       new_mu=(l0*m0[j]+nn[j]*sample_mean.loc[j,:].values)/(l0+nn[j])
#   return new_mu
#
#
#ef get_sigma_map(j, nn, l0, psi0, sample_cov, sample_mean, m0):
#   new_sig=psi0[j]
#   if (j in nn) and nn[j]>2:
#       val=sample_cov.loc[j,:]
#       val=val+l0*nn[j]/(l0+nn[j])*np.outer(sample_mean.loc[j,:]-m0[j], sample_mean.loc[j,:]-m0[j])
#       new_sig=val+psi0[j]
#   return new_sig

        
def percentile(n):
    def percentile_(x):
        return np.percentile(x,n)
    percentile_.__name__= 'percentile_%s' % n
    return percentile_

def load_df(test_id):
    df=pd.read_csv('regeno.diags.7.csv', header=0, index_col=0)
    df1=df[df['var_id']==test_id]
    return df1


#def load_df(var, exclude, sex):
#    test_set = list()
#    for s in var.sample_list:
#        if s in exclude:
#            continue
#        if 'CN' not in var.format_dict:
#            cn=0
#        #if var.info['SVTYPE'] in ['BND']:
#        #    cn=0
#        else:
#            cn = var.genotype(s).get_format('CN')
#            if (var.chrom == 'X' or var.chrom == 'Y') and sex[s] == 1:
#                cn=str(float(cn)*2)
#        test_set.append(vcf_rec(var.var_id, s, var.info['SVTYPE'], var.info['AF'], 
#                    var.genotype(s).get_format('GT'),  cn , var.genotype(s).get_format('AB')))
#
#    test_set = pd.DataFrame(data = test_set, columns=vcf_rec._fields)
#    return test_set


def run_gt_refine(test_id):

#   vcf = Vcf()
#   header = []
#   in_header = True
#   sex={}
#   
#   for line in gender_file:
#       v = line.rstrip().split('\t')
#       sex[v[0]] = int(v[1])
#
#   exclude = []
#   if exclude_file is not None:
#       for line in exclude_file:
#           exclude.append(line.rstrip())
#
#   outf=open(diag_outfile, 'w', 4096)
#   ct=1
#   
#   for line in vcf_in:
#       if in_header:
#           if line[0] == "#":
#              header.append(line)
#              continue
#           else:
#               in_header = False
#               vcf.add_header(header)
#               vcf.add_info('MEDGQR', '1', 'Float', 'Median quality for refined GT')
#               vcf.add_info('Q10GQR', '1', 'Float', 'Q10 quality for refined GT')
#               vcf.add_format('GQR', 1, 'Float', 'Quality of refined genotype.')
#               vcf.add_format('GTR', 1, 'String', 'Refined genotype.')
#               vcf_out.write(vcf.get_header() + '\n')
#
#       v = line.rstrip().split('\t')
#       info = v[7].split(';')
#       svtype = None
#       for x in info:
#           if x.startswith('SVTYPE='):
#               svtype = x.split('=')[1]
#               break
#       # bail if not DEL or DUP prior to reclassification
#       if svtype not in ['DEL', 'BND']:
#           vcf_out.write(line)
#           continue
#       
#       var = Variant(v, vcf)
#       sys.stderr.write("%s\n" % var.var_id)
#       
#       sys.stderr.write("%f\n" % float(var.get_info('AF')))
#       if float(var.get_info('AF'))<0.01:
#           vcf_out.write(line)
#       else:
            df=load_df(test_id)#var, exclude, sex)
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

#   vcf_out.close()
#   vcf_in.close()
#   gender_file.close()
#   outf.close()
#   if exclude_file is not None:
#       exclude_file.close()
              return

def add_arguments_to_parser(parser):
    parser.add_argument('-n', '--id', metavar='<INT>', dest='test_id', type=int, default=None)
 
#
#   parser.add_argument('-i', '--input', metavar='<VCF>', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
#   parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
#   parser.add_argument('-d', '--diag_file', metavar='<STRING>', dest='diag_outfile', type=str, default=None, required=False, help='text file to output method comparisons')
#   parser.add_argument('-g', '--gender', metavar='<FILE>', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
#   parser.add_argument('-e', '--exclude', metavar='<FILE>', dest='exclude', type=argparse.FileType('r'), required=False, default=None, help='list of samples to exclude from classification algorithms')
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'refine genotypes by clustering'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    #with su.InputStream(args.vcf_in) as stream:
    run_gt_refine(args.test_id)

if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
