import argparse, sys, StringIO
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.spatial.distance as ssd
import pysam
sys.path.insert(1,'/gscmnt/gc2802/halllab/abelhj/svtools')
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from collections import namedtuple
import svtools.utils as su
import cProfile , pstats , resource
from scipy.ndimage.interpolation import shift

sys.path.append('/gscmnt/gc2802/halllab/abelhj/svtools/scripts/cncluster_utils')
import CNClusterExact3b_testing_cov2

vcf_rec = namedtuple('vcf_rec', 'varid chr start stop ncarriers sname')

def find_connected_subgraphs(varinfo):

  varinfo.sort_values(['varstart', 'varstop'], ascending=[True, True], inplace=True)
  varinfo.reset_index(inplace=True, drop=True)
  varinfo['comp']=0
  maxpos=0
  clus_id=1
  old_chr=varinfo.chr[0]
  dd=[]
  for ii in range(varinfo.shape[0]):
    if varinfo.chr[ii]!=old_chr:
      sys.stderr.write("Error: merging must be done separately by chromosome\n")
      sys.exit()
    if (maxpos==0) or  (varinfo.varstart[ii]<maxpos):
      maxpos=max(maxpos, varinfo.varstop[ii])
      dd.append(varinfo.index[ii])
    else:
      maxpos=varinfo.varstop[ii]
      varinfo.loc[dd, 'comp']=clus_id
      clus_id=clus_id+1
      dd=[varinfo.index[ii]]
  varinfo.loc[dd, 'comp']=clus_id
  return varinfo

def frac_overlap(a, b):
  overlap=1.0*max(0, min(a[1], b[1]) - max(a[0], b[0]))
  return min(overlap/(a[1]-a[0]), overlap/(b[1]-b[0]))

def min_frac_overlap(a, b, minfrac=0.5):
  ff=frac_overlap(a,b)
  if ff>minfrac:
    return 1
  else:
    return 0
  
def cluster(comp, cn_comp, info_comp, verbose):

  cnag=cn_comp.groupby(['varid']).agg({'cn' : np.var }).reset_index()
  cnag.columns=['varid', 'cnvar']
  print(str(cnag))
  cnag1=cnag.loc[cnag.cnvar>0].copy().reset_index()
  cnvsub2=cn_comp.loc[cn_comp.varid.isin(cnag1.varid)].reset_index(drop=True)
  nvar1=info_comp.loc[info_comp.varid.isin(cnag1.varid)].shape[0]
  cn_comp=0

  if nvar1==0:
    return None
  elif nvar1==1:
    cnvsub2['cluster']=1
    cnvsub2['dist_cluster']=1
  else:
    df2=cnvsub2.pivot_table(index=['varid', 'varstart', 'varstop', 'info_ncarriers'], columns='id', values='cn').reset_index()
    ar=df2.iloc[:, 4:df2.shape[1]].values
    dd=np.linalg.svd(ar, compute_uv=False)
    nvar=1+np.sum(np.cumsum(dd/np.sum(dd))<0.95)
    Z=linkage(ar, method='complete', metric='correlation')
    df2['cluster']=fcluster(Z, nvar, criterion='maxclust')
    df2=df2[['varid', 'varstart', 'varstop', 'info_ncarriers', 'cluster']].copy()
    df2['dist_cluster']=0
    minfrac=0.5
    for clus in df2.cluster.unique():
      temp=df2.loc[df2.cluster==clus]
      if verbose:
        sys.stderr.write(str(clus)+"\n")
        sys.stderr.write("cluster"+str(clus)+"\n")
        print(str(temp))
      if temp.shape[0]==1:
        df2.loc[temp.index, 'dist_cluster']=1
      elif temp.shape[0]==2:
        if min_frac_overlap([temp.varstart[temp.index[0]], temp.varstop[temp.index[0]]], [temp.varstart[temp.index[1]], temp.varstop[temp.index[1]]], minfrac)==1:
          df2.loc[temp.index, 'dist_cluster']=1
        else:
          df2.loc[temp.index, 'dist_cluster']=[1,2]
      else:
        ar=np.empty([temp.shape[0], temp.shape[0]], dtype='float64')
        for ii in range(temp.shape[0]):
          for jj in range(temp.shape[0]):
            ar[ii, jj]=1-frac_overlap([temp.varstart[temp.index[ii]], temp.varstop[temp.index[ii]]], [temp.varstart[temp.index[jj]], temp.varstop[temp.index[jj]]])
        distArray = ssd.squareform(ar)
        Z=linkage(distArray, method='average')
        df2.loc[temp.index, 'dist_cluster']=fcluster(Z, 1-minfrac, criterion='distance')
    cnvsub2=cnvsub2.merge(df2, on=['varid', 'varstart', 'varstop', 'info_ncarriers'])

  gpsub=cnvsub2.drop(['id', 'cn', 'svlen'], axis=1).drop_duplicates()
  gpsub['size']=gpsub['varstop']-gpsub['varstart']
  ag=gpsub.groupby(['comp', 'cluster', 'dist_cluster']).agg({'info_ncarriers' : [np.size,  np.sum ]}).reset_index()
  ag.columns=['comp', 'cluster', 'dist_cluster', 'cluster_nvar', 'cluster_info_ncarriers']
  gpsub=gpsub.merge(ag, on=['comp', 'cluster', 'dist_cluster'])
  gpsub.sort_values(['cluster_info_ncarriers', 'size', 'cluster', 'dist_cluster', 'info_ncarriers'], ascending=[False, True, True, True, False], inplace=True)
  return gpsub


def add_arguments_to_parser(parser):
    parser.add_argument('-l', '--vcf', metavar='<VCF>', dest='lmerged_vcf', help="VCF file containing variants to be output")
    parser.add_argument('-o', '--output', metavar='<STRING>',  dest='outfile', type=str, required=True, help='output file')
    parser.add_argument('-i', '--input', metavar='<STRING>', dest='cnfile', type=str, required=True,  help='copy number file')
    parser.add_argument('-d', '--diag_file', metavar='<STRING>', dest='diag_outfile', type=str,  required=True, help='verbose output file')
    parser.add_argument('-c', '--chrom', metavar='<STRING>', dest='chrom', required=True, help='chrom to analye')
    parser.add_argument('-s', '--samples', metavar='<STRING>', dest='sample_list', required=True, help='list of samples')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')
    parser.add_argument('-r', '--dry_run', dest='dry_run_info_file', metavar='<STRING>', default='',  required=False, help='output info file for dry run')
    parser.add_argument('-t', '--test_comp', dest='test_comp', type=int, required=False, default=-1, help='start at this component for debuggin')

def command_parser():
    parser = argparse.ArgumentParser(description="cross-cohort cnv caller")
    add_arguments_to_parser(parser)
    return parser

def get_info(lmv, chr, sample_list):

  vcf = Vcf()
  in_header = True
  header_lines = list()
  samples=pd.read_csv(sample_list, names=['sampleid'])
  info_set=list()
  nsamps_total=0
  with su.InputStream(lmv) as input_stream:
    for line in input_stream:
      if in_header:
        header_lines.append(line)
        if line[0:9] == '##SAMPLE=':
          nsamps_total+=1
        if line[0:6] == '#CHROM':
          in_header=False
          vcf.add_header(header_lines)
      else:
        v = Variant(line.rstrip().split('\t'), vcf)
        snamestr=v.get_info('SNAME')
        if (v.chrom==chr) and (int(v.info['END'])>v.pos):
          sname=v.info['SNAME']
          # exclude variants with end > start --wtf?
          info_set.append(vcf_rec(v.var_id, v.chrom, v.pos, int(v.info['END']), sname.count(',')+1 , sname))
  info_set = pd.DataFrame(data=info_set, columns=vcf_rec._fields)
  info=info_set.groupby(['chr', 'start', 'stop']).aggregate({'varid':np.min, 'ncarriers':np.sum,  'sname': lambda x : ','.join(x)}).reset_index()

  info.rename(columns={'start':'varstart', 'stop':'varstop', 'ncarriers':'info_ncarriers'}, inplace=True)
  info['svlen']=info['varstop']-info['varstart']
  info_large=info.loc[(info.svlen>10000) & (info.info_ncarriers<20)].copy().reset_index(drop=True)
  info=info.loc[(info.svlen<=10000) | (info.info_ncarriers>=20)].copy().reset_index(drop=True)
  info_large=find_connected_subgraphs(info_large)
  info=find_connected_subgraphs(info)
  info['gp']="info"
  info_large['gp']="info_large"
  info=pd.concat([info, info_large], ignore_index=True)
  info.sort_values(['varstart', 'varstop'], ascending=[True, True], inplace=True)
  info.reset_index(drop=True, inplace=True)
  ord=info[['comp', 'gp']].copy().drop_duplicates()
  ord['rank']=range(ord.shape[0])
  info=info.merge(ord, on=['comp', 'gp'])
  info.drop(['comp', 'gp'], axis=1, inplace=True)
  info.rename(columns={'rank':'comp'}, inplace=True)

  ag_carriers=info[info.info_ncarriers<0.02*nsamps_total].copy().reset_index(drop=True)
  info.drop('sname', axis=1, inplace=True)
  ag_carriers['sname']=ag_carriers['sname'].str.replace(':.', '').str.split(',')
  carriers = pd.DataFrame({
          col:np.repeat(ag_carriers[col].values, ag_carriers['sname'].str.len())
          for col in ag_carriers.columns.drop('sname')}
  ).assign(**{'sname':np.concatenate(ag_carriers['sname'].values)})[ag_carriers.columns]

  carriers=carriers[carriers['sname'].isin(samples['sampleid'])].reset_index(drop=True)
  carriers.rename(columns={'sname':'id'}, inplace=True)
  carriers.drop(['chr', 'varstart', 'varstop', 'svlen'], axis=1, inplace=True)
  return [info, carriers]

def get_cndata(comp, component_pos, info, cntab, verbose):
    temp=component_pos.loc[component_pos.comp==comp].copy().reset_index(drop=True)
    info_comp1=info.loc[info.comp==comp].copy().reset_index(drop=True)
    firstpos=max(0, np.min(temp['varstart'])-10)
    lastpos=min(np.max(temp['varstop'])+10, np.max(component_pos['varstop']))
    if info_comp1.shape[0]<200:
      tabixit=cntab.fetch(component_pos.chr.values[0], firstpos, lastpos)
      if verbose:
        sys.stderr.write("Memory usage info: %s\n" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
      s = StringIO.StringIO("\n".join(tabixit))
      cn=pd.read_table(s,  names=['chr', 'varstart', 'varstop', 'id', 'cn'])
      cn=cn.merge(info_comp1,  on=['chr', 'varstart', 'varstop'])
    else:
       length=lastpos-firstpos+10
       breaks=np.arange(start=firstpos, stop=lastpos+length/5, step=length/10)
       br=pd.DataFrame(np.c_[breaks, shift(breaks, -1)], columns=['start', 'stop'])
       br=br.loc[br.start<lastpos]
       dd=[]
       for ii in range(br.shape[0]):
         if verbose:
           sys.stderr.write("chunk "+str(ii)+"\n")
         tabixit=cntab.fetch(component_pos.chr.values[0], br.start.values[ii], br.stop.values[ii])
         s = StringIO.StringIO("\n".join(tabixit))
         cn1=pd.read_table(s,  names=['chr', 'varstart', 'varstop', 'id', 'cn'])
         cn1=cn1.merge(info_comp1,  on=['chr', 'varstart', 'varstop'])
         if verbose:
           sys.stderr.write("Memory usage info: %s\n" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
         cn1.drop_duplicates(inplace=True)
         dd.append(cn1)
       cn=pd.concat(dd, ignore_index=True)
    cn.drop_duplicates(inplace=True)
    return cn

def prune_info(info_in):
    maxnvar=1000
    component_pos=info_in.groupby(['comp']).agg({'chr': 'first', 'varstart': np.min, 'varstop': np.max, 'varid':np.size}).reset_index()
    cp=component_pos[['comp', 'varid']].copy()
    cp.columns=['comp', 'nvar']
    info1=info_in.merge(cp, on='comp')
    info2=info1.loc[(info1.nvar<maxnvar) | (info1.info_ncarriers>1)].copy()
    cp2=info2.groupby(['comp']).agg({ 'varid':np.size}).reset_index()
    cp2.columns=['comp', 'nvar2']
    info2=info2.merge(cp2, on='comp')
    info3=info2.loc[(info2.nvar2<maxnvar) | (info2.info_ncarriers>2)].copy()
    cp3=info3.groupby(['comp']).agg({ 'varid':np.size}).reset_index()
    cp3.columns=['comp', 'nvar3']
    info3=info3.merge(cp3, on='comp')
    info4=info3.loc[(info3.nvar3<maxnvar) | (info3.info_ncarriers>5)].copy()
    cp4=info4.groupby(['comp']).agg({ 'varid':np.size}).reset_index()
    cp4.columns=['comp', 'nvar4']
    info4=info4.merge(cp4, on='comp')
    info5=info4.loc[(info4.nvar4<maxnvar) | (info4.info_ncarriers>10)].copy()
    cp5=info5.groupby(['comp']).agg({ 'varid':np.size}).reset_index()
    cp5.columns=['comp', 'nvar5']
    info5=info5.merge(cp5, on='comp')
    info=info5[['chr', 'varstart', 'varstop', 'info_ncarriers', 'varid', 'svlen', 'comp']].copy().reset_index(drop=True)
    return info

                                                              
def run_from_args(args):
  cp = cProfile.Profile()
  cp.enable()
  info, carriers=get_info(args.lmerged_vcf, args.chrom, args.sample_list)
  if args.dry_run_info_file!='':
    info.to_csv(args.dry_run_info_file)
    exit(1)
  info=prune_info(info)
  component_pos=info.groupby(['comp']).agg({'chr': 'first', 'varstart': np.min, 'varstop': np.max}).reset_index()

  cntab=pysam.TabixFile(args.cnfile)
  cn=get_cndata(0, component_pos, info, cntab, args.verbose)
  if args.verbose:
    sys.stderr.write("Memory usage info: %s\n" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
  nind=np.unique(cn.loc[cn.comp==0, 'id']).shape[0]
  print(str(nind))
  outf1=open(args.outfile, "w")
  outf2=open(args.diag_outfile, "w")
  header='\t'.join(['#comp', 'cluster', 'dist_cluster', 'start', 'stop', 'nocl', 'bic', 'mean_sep', 'mean_offset', 'cov', 'wts', 'cn_med', 'cn_mad', 'info_ncarriers', 'is_rare', 'mm_corr', 'dist', 'dip_p', 'n_outliers', 'nvar', 'score', 'ptspos', 'ptsend', 'prpos', 'prend', 'is_winner'])
  outf1.write(header+"\n")
  outf2.write(header+"\n")

  for comp in component_pos.comp.unique():
    if(comp==args.test_comp) or (args.test_comp==-1):  
      cn=get_cndata(comp, component_pos, info, cntab, args.verbose)                  
      if args.verbose:
        sys.stderr.write("Memory usage info: %s\n" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
      cn_comp=cn.loc[cn['comp']==comp].copy().reset_index(drop=True)
      info_comp=info.loc[info['comp']==comp].copy().reset_index(drop=True)
      sys.stderr.write("comp="+str(comp)+"\tnvar="+str(info_comp.shape[0])+"\n")
      if (args.verbose):
        print(str(info_comp))
        print( str(info_comp.shape[0]))  
  
      if cn_comp.shape[0]>=nind:
        carriers_comp=carriers.loc[carriers['comp']==comp].copy().reset_index(drop=True)
        region_summary=cluster(comp, cn_comp, info_comp, args.verbose)
        if region_summary is not None:
          for [clus, dist_clus] in region_summary[['cluster', 'dist_cluster']].drop_duplicates().values:
            clus_vars=region_summary.loc[(region_summary.cluster==clus) & (region_summary.dist_cluster==dist_clus)].copy().reset_index(drop=True)
            clus_cn=cn_comp.loc[cn_comp.varid.isin(clus_vars.varid)].copy().reset_index(drop=True)
            clus_carriers=carriers_comp[carriers_comp.varid.isin(clus_vars.varid)].copy().reset_index(drop=True)
            clus=CNClusterExact3b_testing_cov2.CNClusterExact(clus_vars, clus_cn, clus_carriers, args.verbose)
            clus.fit_generic(outf1, outf2)
  
  outf1.close()
  outf2.close()
  cp.disable()
  cp.print_stats()

                                                                                                             
parser=command_parser()
args=parser.parse_args()
print(str(args))
run_from_args(args)

