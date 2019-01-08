import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn import cluster
import gaussian_mixture_constr
import CNWindow1
from statsmodels import robust

class CNClusterExact:
  
  def __init__(self, clus_vars, cn_comp, carriers_comp, nocl_min=1, nocl_max=20, nmads=5):
    self.comp_id=clus_vars.comp[0]
    self.clus_info=clus_vars
    self.cndata=cn_comp
    self.clus_id=clus_vars.cluster[0]
    self.dist_clus_id=clus_vars.dist_cluster[0]
    self.nocl_min=nocl_min
    self.nocl_max=nocl_max
    self.carriers=carriers_comp
    self.nsamp=np.unique(self.cndata.id).size
    self.ncarriers=self.clus_info.cluster_nsamp[0]
    self.nmads=nmads

  def refine_cluster_basic(self, outf2):
    print("refine basic")
    fits=[]
    data=[]
    for ii in range(self.clus_info.shape[0]):
      fit=self.fit_one_window(self.clus_info.varstart.values[ii], self.clus_info.varstop.values[ii],  self.clus_info.nsamp.values[ii])
      fits.append(fit)
    fits=np.concatenate(fits)
    return fits

  def get_chunk_data(self, chunkstart, chunkstop):
    cn=self.cndata[(self.cndata.varstart==chunkstart) & (self.cndata.varstop==chunkstop)].reset_index()
    cn['chunkstart']=chunkstart
    cn['chunkstop']=chunkstop
    cn['clus_id']=self.clus_id
    cn['dist_clus_id']=self.dist_clus_id
    cn.rename(columns={'cn1': 'cn'}, inplace=True)
    return cn

  def fit_one_window(self, chunkstart, chunkstop, nsamp):
    cn1=self.get_chunk_data(chunkstart, chunkstop)
    win=CNWindow1.CNWindow(self.comp_id, self.clus_id, self.dist_clus_id, chunkstart, chunkstop, cn1.cn, self.nocl_min, self.nocl_max)
    fit=win.fit_all_models(nsamp)
    return fit

  def fit_mixture(self, outf2):
    return self.refine_cluster_basic(outf2)

  def fit_generic(self, outf1, outf2):
    #print(str(self.ncarriers)+"\t"+str(self.nsamp))
    #if self.ncarriers>0.0025*self.nsamp:
    if self.ncarriers>=10:
      fit=self.fit_mixture(outf2)
    else:   
      print("try rare")
      fit=self.check_outliers_all(outf2)
    np.savetxt(outf1, fit, delimiter="\t", fmt="%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f")
    return fit

  def check_outliers_all(self, outf2):
    dd=[]
    fits=[]
    for ii in range(self.clus_info.shape[0]):
      cn1=self.get_chunk_data(self.clus_info.varstart.values[ii], self.clus_info.varstop.values[ii])
      cn1['varid']=self.clus_info.varid.values[ii]
      cn1['nsamp_var']=self.clus_info.nsamp.values[ii]
      dd.append(cn1)
    dd=pd.concat(dd)
    cn1=dd.merge(self.carriers,  on=['varid', 'id', 'comp'], how='left')
    cn1['carrier']=np.where(cn1['nsamp'].isnull(), 'non', 'carrier')
    cn1_ag=cn1.groupby(['varid', 'carrier', 'clus_id', 'dist_clus_id', 'comp', 'chunkstart', 'chunkstop', 'nsamp_var'])['cn'].agg({'cn_med': np.median, 'cn_mad':robust.stand_mad}).reset_index()
    print(str(cn1_ag))
    ag_carriers=cn1_ag.loc[cn1_ag.carrier=="carrier"][['varid', 'clus_id', 'dist_clus_id', 'comp', 'cn_med', 'chunkstart', 'chunkstop', 'nsamp_var']].rename(columns={'cn_med':'carrier_med'})
    ag_non=cn1_ag.loc[cn1_ag.carrier=="non"]   
    ag=ag_non.merge(ag_carriers, on=['varid', 'clus_id', 'dist_clus_id', 'comp', 'chunkstart', 'chunkstop', 'nsamp_var'], how='left')
    ag['ll']=ag['cn_med']-self.nmads*ag['cn_mad']
    ag['ul']=ag['cn_med']+self.nmads*ag['cn_mad']
    ag['dist']=(ag['cn_med']-ag['carrier_med'])/ag['cn_mad']
    ag['n_outliers']=0
    ag1=ag.loc[(ag.dist>self.nmads) | (ag.dist < -1*self.nmads) ]
    rare=False
    if ag1.shape[0]>0:
      for varid in ag1['varid'].values:
        cn2=dd.loc[dd.varid == varid]
        cn2=cn2.merge(ag1[['varid', 'll', 'ul', 'cn_med', 'carrier_med', 'cn_mad', 'dist']], on=['varid'])
        if cn2['dist'].values[0]<0:
          cn2=cn2.loc[cn2.cn>cn2.ul]
        else:
          cn2=cn2.loc[cn2.cn<cn2.ll]
        n_outlier=cn2.shape[0]
        ag.loc[ag.varid==varid, 'n_outliers']=n_outlier
        if n_outlier<0.01*self.nsamp and cn2['cn_mad'].values[0]<0.25:
          rare=True
          print("rare variant")
        #res=ag[['comp', 'clus_id',  'chunkstart', 'chunkstop', 'comp', 'carrier_med', 'cn_med', 'cn_mad', 'dist', 'n_outliers', 'nsamp_var', 'nsamp_var', 'nsamp_var', 'nsamp_var']].values
        #res[:,4]=0
        #fits.append(res)
        #if rare:
        #  return res
    #else:
    res=ag[['comp', 'clus_id', 'dist_clus_id', 'chunkstart', 'chunkstop', 'comp', 'carrier_med', 'cn_med', 'cn_mad', 'dist', 'n_outliers', 'nsamp_var', 'nsamp_var', 'nsamp_var', 'nsamp_var']].values
    res[:,5]=0
    if rare:
      return res
    fits.append(res)
    fit=self.fit_mixture(outf2)
    fits.append(fit)
    fits=np.concatenate(fits)
    return fits
    
    
