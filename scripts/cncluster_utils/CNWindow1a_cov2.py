import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.stats import binom, norm, gamma, uniform, dirichlet, shapiro
import gaussian_mixture_constr
from statsmodels import robust
import sys
sys.path.append('/gscmnt/gc2802/halllab/abelhj/svtools/scripts/cncluster_utils')
import dip
#from sklearn import mixture

class CNWindow(object):

  def __init__(self, comp_id, clus_id, clus_dist_id, start, stop, cndata, nocl_min, nocl_max):
    self.data=cndata
    self.nocl_min=nocl_min
    self.nocl_max=nocl_max
    self.clus_id=clus_id
    self.clus_dist_id=clus_dist_id
    self.comp_id=comp_id
    self.start=start
    self.stop=stop
    self.nsamp_rm_outliers=0
    self.nsamp=self.data.size
    self.procdata=np.empty(0)
    self.remove_outliers()
    self.cn_med=np.median(cndata)
    self.cn_mad=robust.stand_mad(cndata)
    
  def remove_outliers(self, qtop=0.0025, qbot=0.9975):
    cn2=np.sort(self.data)
    nsamp=cn2.size
    cn2=cn2+norm.rvs(loc=0, scale=0.005, size=cn2.size)
    self.nsamp_rm_outliers=cn2.size
    self.procdata=cn2.reshape(-1, 1)

  def fit_all_models(self, ncarriers, verbose):
    fits=[]
    bic_col=6
    mm=2*(np.max(self.procdata)-np.min(self.procdata))+2
    if mm>self.nocl_max:
      mm=self.nocl_max
    res=self.fit_one_model(1)
    fits.append(res)
    bic_min=res[bic_col]
    nocl_min=1
    nocl=2
    bic=res[bic_col]
    #ksp=shapiro(self.procdata)[1]
    dipp=dip.diptst1(self.procdata[:,0])[1]
    moveon=False
    if verbose:
      print(str(self.start)+"\t"+str(self.stop)+"\t1\t"+str(bic)+"\t"+str(bic_min)+"\t"+str(dipp))
    while (nocl<mm): # and (bic<=bic_min):
      res=self.fit_one_model(nocl)
      fits.append(res)
      bic=res[bic_col]
      if bic<bic_min:
        bic_min=bic
        nocl_min=nocl
      if verbose:
        print(str(self.start)+"\t"+str(self.stop)+"\t"+str(nocl)+"\t"+str(bic)+"\t"+str(bic_min)+"\t"+str(dipp))
      #save some model fitting
      if nocl_min<=2 and nocl==4:  
        break
      if nocl_min==3 and nocl==6:
        break
      if nocl_min==4 and nocl==8:
        break
      nocl=nocl+1
    if verbose:
      print(str(nocl_min))
    fits=np.vstack(fits)
    if nocl_min>=4:
      fit_temp=pd.DataFrame(fits)
      fit_temp.columns=['comp', 'c1', 'c2', 'start', 'stop', 'nocl', 'bic', 'mm', 'kk', 'vars', 'wts', 'xx', 'yy']
      gmin=np.min(fit_temp.bic)
      min_mm=float(fit_temp.loc[fit_temp.bic==gmin, 'mm'])
      if min_mm<0.52:
        fit_temp['bic1']=-1*(fit_temp['bic']-np.max(fit_temp['bic']))
        a=fit_temp.bic.values
        fit_temp['local_min']=np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
        fit_temp['bic_frac']=fit_temp['bic1']/np.max(fit_temp['bic1'])
        ff=fit_temp.loc[(fit_temp.local_min==True) & (fit_temp.bic_frac>0.7), ['nocl']]
        if ff.shape[0]>0:
          nocl_min=int(np.min(ff['nocl']))
    fits=np.hstack([fits, np.empty([fits.shape[0], 2], dtype='int32'), np.empty([fits.shape[0], 1], dtype='float64')])
    fits[:, 13]=ncarriers
    fits[:,14]=nocl_min
    fits[:,15]=dipp
    return fits

  def fit_one_model(self, nocl):
    gmm=gaussian_mixture_constr.GaussianMixture(n_components=nocl, n_init=10, covariance_type='spherical') 
    gmm.fit(self.procdata)
    cov=','.join( map(str, np.round(gmm.covariances_, 4)))
    wts= ','.join( map(str, np.round(gmm.weights_, 4)))
    bic=gmm.bic(self.procdata)
    aic=gmm.aic(self.procdata)
    lld=gmm.score(self.procdata)
    nn=gmm._n_parameters()
    kk, mm = gmm.get_kk_mm()
    ret=np.array([self.comp_id, self.clus_id, self.clus_dist_id, self.start, self.stop, nocl, bic, mm, kk, cov, wts, self.cn_med, self.cn_mad], dtype='object')
    return ret

    

