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
    #firstel=int(np.floor(nsamp*0.0025))
    #lastel=int(np.floor(nsamp*0.9975))
    #q_002=cn2[firstel]
    #q_998=cn2[lastel]
    #if np.abs(q_002-cn2[0])<0.1:
    #  firstel=0
    #if np.abs(q_998-cn2[nsamp-1])<0.1:
    #  lastel=nsamp-1
    #cn2=cn2[firstel:(lastel+1)]
    cn2=cn2+norm.rvs(loc=0, scale=0.005, size=cn2.size)
    self.nsamp_rm_outliers=cn2.size
    self.procdata=cn2.reshape(-1, 1)

  def fit_all_models(self, ncarriers, verbose):
    fits=[]
    bic_col=6
    mm=3*(np.max(self.procdata)-np.min(self.procdata))
    if mm>self.nocl_max:
      mm=self.nocl_max
    #mm=1+int(2*np.max(self.procdata))
    #if mm<self.nocl_max:
    #  self.nocl_max=mm
    res=self.fit_one_model(1)
    fits.append(res)
    bic_min=res[bic_col]
    cl_min=1
    nocl=2
    bic=res[bic_col]
    ksp=shapiro(self.procdata)[1]
    dipp=dip.diptst(self.procdata[:,0])[1]
    if verbose:
      print(str(self.start)+"\t"+str(self.stop)+"\t1\t"+str(bic)+"\t"+str(bic_min)+"\t"+str(ksp)+"\t"+str(dipp))
    while (nocl<mm): # and (bic<=bic_min):
      res=self.fit_one_model(nocl)
      fits.append(res)
      bic=res[bic_col]
      if verbose:
        print(str(self.start)+"\t"+str(self.stop)+"\t"+str(nocl)+"\t"+str(bic)+"\t"+str(bic_min)+"\t"+str(ksp)+"\t"+str(dipp))
      if bic<bic_min:
        bic_min=bic
      nocl=nocl+1
    cl_min=nocl-1
    if verbose:
      print(str(cl_min))
    fits=np.vstack(fits)
    fits=np.hstack([fits, np.empty([fits.shape[0], 2], dtype='int32')])
    fits[:, 13]=ncarriers
    fits[:,14]=cl_min
    return fits

  def fit_one_model(self, nocl):
    gmm=gaussian_mixture_constr.GaussianMixture(n_components=nocl, n_init=10, covariance_type='spherical') 
    gmm.fit(self.procdata)
    cov=','.join( map(str, np.round(gmm.covariances_, 4)))
    wts= ','.join( map(str, np.round(gmm.weights_, 4)))
    #print(str(self.start)+"\t"+str(self.stop)+"\t"+str(nocl))
    #print(str(gmm.covariances_))
    #print(str(gmm.weights_))
    #print(str(gmm.means_))
    #cov=np.array_str(gmm.covariances_)
    #wts=np.array_str(gmm.weights_)
    bic=gmm.bic(self.procdata)
    aic=gmm.aic(self.procdata)
    lld=gmm.score(self.procdata)
    nn=gmm._n_parameters()
    #print(str(self.start)+"\t"+str(self.stop)+"\t"+str(nocl)+"\t"+str(bic)+"\t"+str(lld)+"\t"+str(nn)+"\t"+str(gmm.means_)+"\t"+str(gmm.covariances_)+"\t"+str(gmm.weights_)) 
    kk, mm = gmm.get_kk_mm()
    ret=np.array([self.comp_id, self.clus_id, self.clus_dist_id, self.start, self.stop, nocl, bic, mm, kk, cov, wts, self.cn_med, self.cn_mad], dtype='object')
    return ret

    

