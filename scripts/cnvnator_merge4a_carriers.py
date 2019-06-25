import argparse, sys, StringIO
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.spatial.distance as ssd
import pysam
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from collections import namedtuple
import svtools.utils as su
import cProfile , pstats

import CNClusterExact3a

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
      print("Error: merging must be done separately by chromosome\n")
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

  cnvsub1=cn_comp.merge(info_comp, on=['comp', 'chr', 'varstart', 'varstop'])
  cnvsub1.rename(columns={'cn1':'cn'}, inplace=True)
  cnag=cnvsub1.groupby(['varid']).agg({'cn' : np.var }).reset_index()
  cnag.columns=['varid', 'cnvar']
  cnag1=cnag.loc[cnag.cnvar>0].copy().reset_index()
  cnvsub2=cnvsub1.merge(cnag1, on=['varid']).reset_index()
  info_comp=info_comp.merge(cnag1, on=['varid'])
  
  if info_comp.shape[0]==1:
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
        print(str(clus))
        print("cluster"+str(clus)+"\n")
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
  return cnvsub2[['varid', 'id', 'cn', 'varstart', 'varstop','comp', 'info_ncarriers', 'cluster', 'dist_cluster']].copy()


def summarize_clusters(cn_clustered):  
  gpsub=cn_clustered.drop(['id', 'cn'], axis=1).drop_duplicates()
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
        if v.chrom==chr:
          sname=v.info['SNAME']
          info_set.append(vcf_rec(v.var_id, v.chrom, v.pos, int(v.info['END']), sname.count(',')+1 , sname))
  info_set = pd.DataFrame(data=info_set, columns=vcf_rec._fields)
  info=info_set.groupby(['chr', 'start', 'stop']).aggregate({'varid':np.min, 'ncarriers':np.sum,  'sname': lambda x : ','.join(x)}).reset_index()

  info.rename(columns={'start':'varstart', 'stop':'varstop', 'ncarriers':'info_ncarriers'}, inplace=True)
  info['svlen']=info['varstop']-info['varstart']
  info_large=info.loc[(info.svlen>100000) & (info.info_ncarriers<20)].copy().reset_index(drop=True)
  info=info.loc[(info.svlen<=100000) | (info.info_ncarriers>=20)].copy().reset_index(drop=True)
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
                                                                                                             

cp = cProfile.Profile()
cp.enable()

parser=command_parser()
args=parser.parse_args()

info, carriers=get_info(args.lmerged_vcf, args.chrom, args.sample_list)
component_pos=info.groupby(['comp']).agg({'chr': 'first', 'varstart': np.min, 'varstop': np.max}).reset_index()


cntab=pysam.TabixFile(args.cnfile)
if component_pos.shape[0]<100:
  tabixit=cntab.fetch( component_pos.chr.values[0], max(0, component_pos.varstart.values[0]-10000), np.max(component_pos.varstop))
else :
  tabixit=cntab.fetch( component_pos.chr.values[0], max(0, component_pos.varstart.values[0]-10000), min(component_pos.varstop.values[100]+10000, np.max(component_pos.varstop)))
  
s = StringIO.StringIO("\n".join(tabixit))
cn=pd.read_table(s,  names=['chr', 'varstart', 'varstop', 'id', 'cn1'])
cn=cn.merge(info,  on=['chr', 'varstart', 'varstop'])
cn.drop([ 'varid', 'info_ncarriers'], axis=1, inplace=True)
cn.drop_duplicates(inplace=True)
nind=np.unique(cn.loc[cn.comp==0, 'id']).shape[0]


outf1=open(args.outfile, "w")
outf2=open(args.diag_outfile, "w")


for comp in component_pos.comp.unique():
  if (comp>0) and (comp%100==0):
    lastcomp=min(comp+100, component_pos.shape[0]-1)
    tabixit=cntab.fetch( component_pos.chr.values[comp], max(0, component_pos.varstart.values[comp]-10000), min(component_pos.varstop.values[lastcomp]+10000, np.max(component_pos.varstop)))
    s = StringIO.StringIO("\n".join(tabixit))
    cn=pd.read_table(s,  names=['chr', 'varstart', 'varstop', 'id', 'cn1'])
    cn=cn.merge(info,  on=['chr', 'varstart', 'varstop'])
    cn.drop([ 'varid', 'info_ncarriers'], axis=1, inplace=True)
    cn.drop_duplicates(inplace=True)
                      
  print("comp="+str(comp))
  cn_comp=cn.loc[cn['comp']==comp].copy().reset_index(drop=True)
  info_comp=info.loc[info['comp']==comp].copy().reset_index(drop=True)
  if (args.verbose):
    print(str(info_comp))
    print( str(info_comp.shape[0]))
  
  
  if cn_comp.shape[0]>=nind:
    if info_comp.shape[0]>10000:
      info_comp=info_comp.loc[info_comp.info_ncarriers>10].copy().reset_index(drop=True)
    elif info_comp.shape[0]>2000:
      info_comp=info_comp.loc[info_comp.info_ncarriers>5].copy().reset_index(drop=True)
    elif info_comp.shape[0]>500:
      info_comp=info_comp.loc[info_comp.info_ncarriers>1].copy().reset_index(drop=True)
    carriers_comp=carriers.loc[carriers['comp']==comp].copy().reset_index(drop=True)
    cn_region=cluster(comp, cn_comp, info_comp, args.verbose)
    region_summary=summarize_clusters(cn_region)
    for [clus, dist_clus] in region_summary[['cluster', 'dist_cluster']].drop_duplicates().values:
      clus_vars=region_summary.loc[(region_summary.cluster==clus) & (region_summary.dist_cluster==dist_clus)].copy().reset_index(drop=True)
      clus=CNClusterExact3a.CNClusterExact(clus_vars, cn_comp, carriers_comp, args.verbose)
      clus.fit_generic(outf1, outf2)
  
outf1.close()
outf2.close()

cp.disable()
cp.print_stats()
