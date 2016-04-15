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


CN_rec = namedtuple ('CN_rec', 'var_id sample svtype svlen AF GT CN AB log_len log2r')



# http://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


def to_bnd_strings(var, fixed_gts):

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
    var1=var.get_var_string(fixed_gts)

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
    var2=var.get_var_string(fixed_gts)
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
            curr_rec[1] = next_rec[1]
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
            frac_overlap = reciprocal_overlap([var_start, var_end], class_overlap[me_class])
            if frac_overlap > best_frac_overlap:
                best_frac_overlap = frac_overlap
                best_feature = me_class


        if best_frac_overlap >= threshold:
            return best_feature

    return None

def lowQuantile(xx):
    return np.percentile(xx,2.5)

def highQuantile(xx):
    return np.percentile(xx,97.5)

def lld(xx, mean, sd):
    ll = 1 / sd * math.exp(-(xx-mean) * (xx-mean) / (2*sd*sd))
    return ll

def calc_params(vcf_path):

    tSet = list()
    epsilon=0.1
    header=[]
    
    in_header = True
    vcf = Vcf()
    if vcf_path.endswith('.gz'):
        vcf_file = gzip.open(vcf_path, 'rb')
    else:
        vcf_file = open(vcf_path, 'r')

    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line)
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                    in_header = False
                    vcf.add_header(header)
                continue
        else:
            v = line.rstrip().split('\t')
            info = v[7].split(';')
            svtype = None
            for x in info:
                if x.startswith('SVTYPE='):
                    svtype = x.split('=')[1]
                    break

            if svtype not in ['DEL', 'DUP'] or v[0]=="X" or v[0]=="Y":
                continue

            var = Variant(v, vcf)
    
            for sample in vcf_samples:
                sample_genotype = var.genotype(sample)
                if sample_genotype.get_format('GT') != './.':
                    log2r = math.log((float(sample_genotype.get_format('CN'))+ epsilon)/2,2)  #to avoid log(0)
                    tSet.append(CN_rec(var.var_id, sample, var.info['SVTYPE'], abs(float(var.info['SVLEN'])), var.info['AF'],
                        sample_genotype.get_format('GT'),  sample_genotype.get_format('CN'), sample_genotype.get_format('AB'), math.log(abs(float(var.info['SVLEN']))), log2r))

    df=pd.DataFrame(tSet, columns=CN_rec._fields)
    #exclude from training data, DELs and DUPs with CN in the tails of the distribution
    df['q_low']=df.groupby(['sample', 'svtype', 'GT'])['log2r'].transform(lowQuantile)
    df['q_high']=df.groupby(['sample', 'svtype', 'GT'])['log2r'].transform(highQuantile)
    df=df[(df.log2r>=df.q_low) & (df.log2r<=df.q_high)]
    #df.to_csv('./train.csv')

    #adjust copy number for small deletions (<1kb), no strong relationship b/w cn and size for dups evident so far
    small_het_dels = df[(df.svtype=="DEL") & (df.GT=="0/1") & (df.svlen<1000) & (df.svlen>=50)]
    small_hom_dels = df[(df.svtype=="DEL") & (df.GT=="1/1") & (df.svlen<1000) & (df.svlen>=50)]
    het_del_mean=np.mean(df[(df.svlen>1000) & (df.GT=="0/1") & (df.svtype=="DEL")]['log2r'])
    hom_del_mean=np.mean(df[(df.svlen>1000) & (df.GT=="1/1") & (df.svtype=="DEL")]['log2r'])
    small_het_dels['offset']=small_het_dels['log2r']-het_del_mean
    small_hom_dels['offset']=small_hom_dels['log2r']-hom_del_mean
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        hom_del_fit=smf.ols('offset~log_len',small_hom_dels).fit()
        het_del_fit=smf.ols('offset~log_len',small_het_dels).fit()
        #print hom_del_fit.summary()
        #print het_del_fit.summary()
        small_hom_dels['log2r_adj'] = small_hom_dels['log2r'] - hom_del_fit.predict(small_hom_dels)
        small_het_dels['log2r_adj'] = small_het_dels['log2r'] - het_del_fit.predict(small_het_dels)

    small_dels=small_hom_dels.append(small_het_dels)
    small_dels=small_dels[['var_id', 'sample', 'svtype', 'svlen', 'AF', 'GT', 'CN', 'log_len', 'log2r', 'q_low', 'q_high', 'log2r_adj']]

    # dels of length<100 bp are excluded here
    df1=df[(df.svtype!="DEL") | (df.GT=="0/0") | (df.svlen>=1000)]
    df1['log2r_adj']=df1['log2r']
    df1=df1.append(small_dels)

    params=df1.groupby(['sample', 'svtype', 'GT'])['log2r_adj'].aggregate([np.mean,np.var, len]).reset_index()
    params=pd.pivot_table(params, index=['sample', 'svtype'], columns='GT', values=['mean', 'var', 'len']).reset_index()    
    params.columns=['sample', 'svtype', 'mean0', 'mean1', 'mean2', 'var0', 'var1', 'var2', 'len0', 'len1', 'len2']
    params['std_pooled']=np.sqrt((params['var0']*params['len0']+params['var1']*params['len1']+params['var2']*params['len2'])/(params['len0']+params['len1']+params['len2']))
    #params.to_csv('./params.csv')
    return (params, het_del_fit, hom_del_fit)


def rd_support_nb(temp, p_cnv):
    tr = pd.DataFrame({'p0' : [1.0, 0.1, 0.0], 'p1' : [0.0, 0.7, 0.25], 'p2' : [0.0, 0.2, 0.75], 'GT' : ["0/0", "0/1", "1/1"]})
    temp = pd.merge(temp, tr, on='GT', how='left')
    temp['p_mix'] = temp['lld0'] * temp['p0'] + temp['lld1'] * temp['p1'] + temp['lld2'] * temp['p2']
    return np.log(p_cnv)+np.sum(np.log(temp['p_mix'])) > np.log(1-p_cnv)+np.sum(np.log(temp['lld0']))
   

def has_rd_support_by_nb(test_set, het_del_fit, hom_del_fit, params, p_cnv = 0.5):
    svtype=test_set['svtype'][0]
    svlen=test_set['svlen'][0]
    log_len=test_set['log_len'][0]
    
    if svtype == 'DEL' and svlen<1000:
        params1=params[params.svtype=='DEL'].copy()
        if svlen<50:
            params1['log_len']=math.log(50)
        else:
            params1['log_len']=log_len
        params1['mean1_adj'] = params1['mean1'] + het_del_fit.predict(params1)
        params1['mean2_adj'] = params1['mean2'] + hom_del_fit.predict(params1)
    else:
        params1=params.copy()
        params1['mean1_adj'] = params1['mean1']
        params1['mean2_adj'] = params1['mean2']

    v0=test_set[test_set.GT=="0/0"]['log2r'].values
    v1=test_set[test_set.GT=="0/1"]['log2r'].values
    v2=test_set[test_set.GT=="1/1"]['log2r'].values

    if len(v0)>0:
        med0=np.median(v0)
    else:
        if len(v1)>0:
            med0=med1=np.median(v1)
        elif len(v2)>0:
            med0=med1=med2=np.median(v2)
        else:
            return False

    if len(v1)>0:
        med1=np.median(v1)
    else:
        med1=med0
    if len(v2)>0:
        med2=np.median(v2)
    else:
        med2=med1

    if svtype=='DEL' and ( med1>med0 or med2>med0 ):
        return False
    elif svtype=='DUP' and (med1<med0 or med2<med0):
        return False

    mm=pd.merge(test_set, params1, how='left')

    mm['lld0'] = mm.apply(lambda row:lld(row["log2r"], row["mean0"],row["std_pooled"]), axis=1)
    mm['lld1'] = mm.apply(lambda row:lld(row["log2r"], row["mean1_adj"],row["std_pooled"]), axis=1)
    mm['lld2'] = mm.apply(lambda row:lld(row["log2r"], row["mean2_adj"],row["std_pooled"]), axis=1)
   
    return  rd_support_nb(mm, p_cnv)


def load_df(var, exclude, sex):
    
    epsilon=0.1
    test_set = list()

    for s in var.sample_list:
        if s in exclude:
            continue
        cn = var.genotype(s).get_format('CN')
        if (var.chrom == 'X' or var.chrom == 'Y') and sex[s] == 1:
            cn=str(float(cn)*2)
        log2r = math.log((float(cn)+epsilon)/2, 2)  # to avoid log(0)
        test_set.append(CN_rec(var.var_id, s, var.info['SVTYPE'], abs(float(var.info['SVLEN'])), var.info['AF'],
             var.genotype(s).get_format('GT'),  cn , var.genotype(s).get_format('AB'), math.log(abs(float(var.info['SVLEN']))), log2r))

    test_set = pd.DataFrame(data = test_set, columns=CN_rec._fields)
    return test_set

# test for read depth support of low frequency variants
def has_low_freq_depth_support(test_set, mad_threshold=2, absolute_cn_diff=0.5):

    mad_quorum = 0.5 # this fraction of the pos. genotyped results must meet the mad_threshold
    
    hom_ref_cn=test_set[test_set.GT=="0/0"]['CN'].values.astype(float)
    hom_het_alt_cn=test_set[(test_set.GT=="0/1") | (test_set.GT=="1/1")]['CN'].values.astype(float)

    if len(hom_ref_cn) > 0:
        cn_median = np.median(hom_ref_cn)
        cn_mad = mad(hom_ref_cn)
    else:
        cn_median = None
        cn_mad = None

    # bail after writing out diagnostic info, if no ref samples or all ref samples
    if (len(hom_ref_cn) == 0 or len(hom_het_alt_cn) == 0):
        return False

    # tally up the pos. genotyped samples meeting the mad_threshold

    resid=hom_het_alt_cn-cn_median
    if test_set['svtype'][0]=='DEL':
        resid=-resid
    
    resid=resid[(resid > (cn_mad * mad_threshold) ) & (resid>absolute_cn_diff)]

    if float(len(resid))/len(hom_het_alt_cn)>mad_quorum:
        return True
    else:
        return False

# test whether variant has read depth support by regression
def has_high_freq_depth_support(df, slope_threshold, rsquared_threshold):
    
    rd = df[[ 'AB', 'CN']][df['AB']!='.'].values.astype(float)
    if len(np.unique(rd[0,:])) > 1 and len(np.unique(rd[1,:])) > 1:
        
        (slope, intercept, r_value, p_value, std_err) = stats.linregress(rd)
        if df['svtype'][0] == 'DEL':
            slope=-slope
        #sys.stderr.write(df['var_id'][0]+"\t"+str(slope)+"\t"+str(r_value)+"\n")

        if (slope < slope_threshold or r_value*r_value < rsquared_threshold):
            return False
        return True
    return False

def has_rd_support_by_ls(df, slope_threshold, rsquared_threshold, num_pos_samps, mad_threshold=2, absolute_cn_diff=0.5):

    min_pos_samps_for_regression=10
    if num_pos_samps>min_pos_samps_for_regression:
        return has_high_freq_depth_support(df, slope_threshold, rsquared_threshold)
    else:
        return has_low_freq_depth_support(df, mad_threshold, absolute_cn_diff)
    return False

def has_rd_support_hybrid(df, het_del_fit, hom_del_fit, params, p_cnv, slope_threshold, rsquared_threshold, num_pos_samps):
  
    hybrid_support=False
    nb_support=has_rd_support_by_nb(df, het_del_fit, hom_del_fit, params, p_cnv)
    ls_support=has_rd_support_by_ls(df, slope_threshold, rsquared_threshold, num_pos_samps)

    if nb_support and ls_support:
        hybrid_support=True
    elif nb_support and has_rd_support_by_ls(df, 2*slope_threshold, 2*rsquared_threshold, num_pos_samps, 2, 0.75):
        hybrid_support=True
    elif ls_support and has_rd_support_by_nb(df, het_del_fit, hom_del_fit, params, 0.2*p_cnv):
        hybrid_support=True
    return [ls_support, nb_support, hybrid_support]


# primary function
def sv_classify(vcf_in, vcf_out, gender_file, exclude_file, ae_dict, f_overlap, slope_threshold, rsquared_threshold, p_cnv, het_del_fit, hom_del_fit, params, diag_outfile, method):

    vcf = Vcf()
    header = []
    in_header = True
    sex = {}
    # read sample genders
    for line in gender_file:
        v = line.rstrip().split('\t')
        sex[v[0]] = int(v[1])

    exclude = []
    if exclude_file is not None:
        for line in exclude_file:
            exclude.append(line.rstrip())

    if diag_outfile is not None:
        outf=open(diag_outfile, 'w', 4096)
        outf.write("varid\torig_svtype\tsvlen\tnum_pos_samps\tnb_support\tls_support\thybrid_support\thas_rd_support\n")

    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf_out.write(vcf.get_header() + '\n')

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
        
        var = Variant(v, vcf)

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

        #count positively genotyped samples
        num_pos_samps = 0
        num_total_samps=len(var.sample_list)

        for s in var.sample_list:
            if var.genotype(s).get_format('GT') not in ["./.", "0/0"]:
                num_pos_samps += 1

        nb_support = False
        ls_support = False
        hybrid_support = False
        has_rd_support = False

        if num_pos_samps == 0:
            vcf_out.write(line)
        else:
            df=load_df(var, exclude, sex)
            if method=='large_sample':
                ls_support = has_rd_support_by_ls(df, slope_threshold, rsquared_threshold, num_pos_samps)
                has_rd_support=ls_support
            elif method=='naive_bayes':
                nb_support = has_rd_support_by_nb(df, het_del_fit, hom_del_fit, params, p_cnv)
                has_rd_support=nb_support
            elif method=='hybrid':
                ls_support, nb_support, hybrid_support = has_rd_support_hybrid(df, het_del_fit, hom_del_fit, params, p_cnv, slope_threshold, rsquared_threshold, num_pos_samps)
                has_rd_support=hybrid_support

            if has_rd_support:
               vcf_out.write(line)
            else:
                for m_var in to_bnd_strings(var, True):
                    vcf_out.write(m_var + '\n')

            if diag_outfile is not None:
              svlen=df['svlen'][0]
              outf.write(var.var_id+"\t"+svtype+"\t"+str(svlen)+"\t"+str(num_pos_samps)+"\t"+str(nb_support)+"\t"+str(ls_support)+"\t"+str(hybrid_support)+"\t"+str(has_rd_support)+"\n")

    vcf_out.close()
    if diag_outfile is not None:
        outf.close()
    vcf_in.close()
    vcf_out.close()
    gender_file.close()
    if exclude_file is not None:
        exclude_file.close()

    return


def get_ae_dict(ae_path):
    if ae_path.endswith('.gz'):
        ae_bedfile = gzip.open(ae_path, 'rb')
    else:
        ae_bedfile = open(ae_path, 'r')
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
    ae_bedfile.close()
    return ae_dict

def run_reclassifier(vcf_file, vcf_out, sex_file, ae_path, f_overlap, exclude_list, slope_threshold, rsquared_threshold, training_data, method, diag_outfile):

    ae_dict = None
    params = None
    het_del_fit = None
    hom_del_fit = None
    p_cnv=0.5       # prior probability that CNV is real

    if ae_path is not None:
        sys.stderr.write("loading annotations\n")
        ae_dict=get_ae_dict(ae_path)
    
    if(method!="large_sample"):
          sys.stderr.write("calculating parameters\n")
          #calculate per-sample CN profiles on training set
          [params, het_del_fit, hom_del_fit]=calc_params(training_data)

    sys.stderr.write("reclassifying\n")
    sv_classify(vcf_file,
                vcf_out,
                sex_file,
                exclude_list,
                ae_dict,
                f_overlap,
                slope_threshold,
                rsquared_threshold,
                p_cnv,
                het_del_fit,
                hom_del_fit,
                params, 
                diag_outfile,
                method)

def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--input', metavar='<VCF>', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-o', '--output', metavar='<VCF>', dest='vcf_out', type=argparse.FileType('w'), default=sys.stdout, help='VCF output [stdout]')
    parser.add_argument('-g', '--gender', metavar='<FILE>', dest='gender', type=argparse.FileType('r'), required=True, default=None, help='tab delimited file of sample genders (male=1, female=2)\nex: SAMPLE_A\t2')
    parser.add_argument('-a', '--annotation', metavar='<BED>', dest='ae_path', type=str, default=None, help='BED file of annotated elements')
    parser.add_argument('-f', '--fraction', metavar='<FLOAT>', dest='f_overlap', type=float, default=0.9, help='fraction of reciprocal overlap to apply annotation to variant [0.9]')
    parser.add_argument('-e', '--exclude', metavar='<FILE>', dest='exclude', type=argparse.FileType('r'), required=False, default=None, help='list of samples to exclude from classification algorithms')
    parser.add_argument('-s', '--slope_threshold', metavar='<FLOAT>', dest='slope_threshold', type=float, default=1.0, help='minimum slope absolute value of regression line to classify as DEL or DUP[1.0]')
    parser.add_argument('-r', '--rsquared_threshold', metavar='<FLOAT>', dest='rsquared_threshold', type=float, default=0.2, help='minimum R^2 correlation value of regression line to classify as DEL or DUP [0.2], for large sample reclassification')
    parser.add_argument('-t', '--tSet', metavar='<STRING>', dest='tSet', type=str, default=None, required=False, help='high quality deletions & duplications training dataset[vcf], required by naive Bayes reclassification')
    parser.add_argument('-m', '--method', metavar='<STRING>', dest='method', type=str, default="large_sample", required=False, help='reclassification method, one of (large_sample, naive_bayes, hybrid)', choices=['large_sample', 'naive_bayes', 'hybrid'])
    parser.add_argument('-d', '--diag_file', metavar='<STRING>', dest='diag_outfile', type=str, default=None, required=False, help='text file to output method comparisons')
    parser.set_defaults(entry_point=run_from_args)

def description():
    return 'reclassify DEL and DUP based on read depth information'

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    if args.tSet is None:
        if args.method!="large_sample":
            sys.stderr.write("Training data required for naive Bayes or hybrid classifiers\n")
            parser.print_help()
            sys.exit(1)
    with su.InputStream(args.vcf_in) as stream:
        run_reclassifier(stream, args.vcf_out, args.gender, args.ae_path, args.f_overlap, args.exclude, args.slope_threshold, args.rsquared_threshold, args.tSet, args.method, args.diag_outfile)

if __name__ == '__main__':
    parser = command_parser()
    args=parser.parse_args()
    sys.exit(args.entry_point(args))
