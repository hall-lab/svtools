import argparse, sys, StringIO
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.spatial.distance as ssd

sys.path.append('/gscmnt/gc2802/halllab/sv_aggregate/dev/svtools/scripts/cncluster_utils')
import CNClusterExact1

def find_connected_subgraphs(varinfo):
  varinfo.reset_index(drop=True, inplace=True)
  varinfo=varinfo.sort_values(['varstart', 'varstop'], ascending=[True, True]).reset_index()
  varinfo=varinfo.groupby(['chr', 'varstart', 'varstop']).agg({'varid': np.min, 'nsamp': np.sum}).reset_index()
  varinfo['comp']=0
  maxpos=0
  clus_id=1
  dd=[]
  for ii in range(varinfo.shape[0]):
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
  
def cluster(comp, cn_comp, info_comp, carriers_comp):
  cnvsub=cn_comp.merge(info_comp, on=['comp', 'chr', 'varstart', 'varstop'])
  cnvsub1=cnvsub.merge(carriers_comp, on=['varid', 'id', 'comp'], how='left')
  cnvsub1['carrier']=np.where(cnvsub1['nsamp_y'].isnull(), 'non', 'carrier')
  cnvsub1.rename(columns={'nsamp_x':'nsamp', 'cn1':'cn'}, inplace=True)
  if info_comp.shape[0]==1:
    cnvsub1['cluster']=1
    cnvsub1['dist_cluster']=1
  else:
    df2=cnvsub1.pivot_table(index=['varid', 'varstart', 'varstop', 'nsamp'], columns='id', values='cn').reset_index()
    ar=df2.iloc[:, 5:df2.shape[1]].values
    dd=np.linalg.svd(ar, compute_uv=False)
    nvar=1+np.sum(np.cumsum(dd/np.sum(dd))<0.95)
    Z=linkage(ar, method='complete', metric='correlation')
    df2['cluster']=fcluster(Z, nvar, criterion='maxclust')
    df2=df2[['varid', 'varstart', 'varstop', 'nsamp', 'cluster']]
    df2['dist_cluster']=0
    minfrac=0.5
    for clus in df2.cluster.unique():
      print(str(clus))
      temp=df2.loc[df2.cluster==clus]
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
    cnvsub1=cnvsub1.merge(df2, on=['varid', 'varstart', 'varstop', 'nsamp'])
  return cnvsub1[['varid', 'id', 'cn', 'varstart', 'varstop','comp', 'nsamp', 'carrier', 'cluster', 'dist_cluster']]


def summarize_clusters(cn_clustered):
  gpsub=cn_clustered[['varid', 'nsamp', 'varstart', 'varstop', 'comp', 'cluster', 'dist_cluster']].drop_duplicates()
  gpsub['size']=gpsub['varstop']-gpsub['varstart']
  ag=gpsub.groupby(['comp', 'cluster', 'dist_cluster']).agg({'nsamp' : [np.size,  np.sum ] }).reset_index()
  ag.columns=['comp', 'cluster', 'dist_cluster', 'num_var', 'cluster_nsamp']
  gpsub=gpsub.merge(ag, on=['comp', 'cluster', 'dist_cluster'])
  rare_clus=gpsub.loc[gpsub.cluster_nsamp<0.005*nind]
  gpsub.loc[rare_clus.index, 'dist_cluster']=0
  gpsub=gpsub[['varid', 'nsamp', 'varstart', 'varstop', 'comp', 'cluster', 'dist_cluster',  'size']]
  ag=gpsub.groupby(['comp', 'cluster', 'dist_cluster']).agg({'nsamp' : [np.size,  np.sum ] }).reset_index()
  ag.columns=['comp', 'cluster', 'dist_cluster', 'num_var', 'cluster_nsamp']
  gpsub=gpsub.merge(ag, on=['comp', 'cluster', 'dist_cluster'])
  gpsub['size']=gpsub['size']*gpsub['nsamp']/gpsub['cluster_nsamp']
  gp=gpsub.groupby(['cluster', 'dist_cluster']).agg({'size': np.sum, 'varstart':np.min, 'varstop':np.max}).reset_index()
  gp.rename(columns={'varstart':'cluster_start', 'varstop':'cluster_stop', 'size':'cluster_mean_size'}, inplace=True)
  gpsub=gpsub.merge(gp, on=['cluster', 'dist_cluster'])
  gpsub.sort_values(['cluster_nsamp', 'cluster_mean_size', 'cluster', 'dist_cluster', 'nsamp'], ascending=[False, True, True, True, False], inplace=True)
  return gpsub

def add_arguments_to_parser(parser):
    parser.add_argument('-c', '--carriers', metavar='<STRING>', dest='carriers_file', type=str, 
default='/Users/habel/halllab/sv_aggregate/cnvnator/vs_gs/carriers.rare.txt.gz', help='carriers file')
    parser.add_argument('-o', '--output', metavar='<STRING>',  dest='outfile', type=str, default='xx.txt', help='output file')
    parser.add_argument('-i', '--input', metavar='<STRING>', dest='cnfile', type=str, 
default='/Users/habel/halllab/sv_aggregate/cnvnator/vs_gs/cnvnator.all.1.chr6.txt.gz', help='copy number file')
    parser.add_argument('-f', '--info', metavar='<STRING>', dest='info_file', type=str, 
default='/Users/habel/halllab/sv_aggregate/cnvnator/vs_gs/clusters.chr6.txt', help='variant info file')
    parser.add_argument('-d', '--diag_file', metavar='<STRING>', dest='diag_outfile', type=str, default='yy.txt', required=False, help='verbose output file')

def command_parser():
    parser = argparse.ArgumentParser(description="cross-cohort cnv caller")
    add_arguments_to_parser(parser)
    return parser

parser=command_parser()
args=parser.parse_args()

info=pd.read_table(args.info_file, names=['chr', 'varstart', 'varstop', 'varid', 'nsamp'])
info=info.groupby(['chr', 'varstart', 'varstop']).agg({'varid': np.min, 'nsamp': np.sum}).reset_index()
info['svlen']=info['varstop']-info['varstart']
info_large=info.loc[(info.svlen>100000) & (info.nsamp<20)]
info=info.loc[(info.svlen<=100000) | (info.nsamp>=20)]
info_large=find_connected_subgraphs(info_large)
info=find_connected_subgraphs(info)
info['gp']="info"
info_large['gp']="info_large"
info=pd.concat([info, info_large], ignore_index=True)
info=info.sort_values(['varstart', 'varstop'], ascending=[True, True])
info.reset_index(drop=True, inplace=True)
ord=info[['comp', 'gp']].drop_duplicates()
ord['rank']=range(ord.shape[0])
info=info.merge(ord, on=['comp', 'gp'])
info=info[['chr', 'varstart', 'varstop', 'varid', 'nsamp', 'rank']].rename(columns={'rank':'comp'})

carriers=pd.read_table(args.carriers_file, names=['chr', 'varstart', 'varstop', 'id'])
cn=pd.read_table(args.cnfile, names=["id", "cn1", "cn2", "chr", "varstart", "varstop", "svlen", "varid", "nsamp", "comp"])[['id', 'cn1', 'cn2', 'chr', 'varstart', 'varstop']]                       

carriers=carriers.merge(info, on=['chr', 'varstart', 'varstop'])[['varid', 'id', 'nsamp', 'comp']].drop_duplicates()
cn=cn.merge(info, on=['chr', 'varstart', 'varstop'])[['id', 'cn1', 'chr', 'varstart', 'varstop', 'comp']].drop_duplicates()
nind=cn.loc[cn.comp==0]['id'].unique().size

outf1=open(args.outfile, "w")
outf2=open(args.diag_outfile, "w")


for comp in cn.comp.unique():
  print("comp="+str(comp))
  cn_comp=cn.loc[cn['comp']==comp]
  info_comp=info.loc[info['comp']==comp]
  print( str(info_comp.shape[0]))
  if info_comp.shape[0]>500:
    info_comp=info_comp.loc[info_comp.nsamp>1]
  carriers_comp=carriers.loc[carriers['comp']==comp]
  cn_region=cluster(comp, cn_comp, info_comp, carriers_comp)
  region_summary=summarize_clusters(cn_region)
  for [clus, dist_clus] in region_summary[['cluster', 'dist_cluster']].drop_duplicates().values:
    clus_vars=region_summary.loc[(region_summary.cluster==clus) & (region_summary.dist_cluster==dist_clus)].reset_index()
    clus=CNClusterExact1.CNClusterExact(clus_vars, cn_comp, carriers_comp)
    clus.fit_generic(outf1, outf2)
  

outf1.close()
outf2.close()
