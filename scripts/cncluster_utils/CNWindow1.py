import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.stats import binom, norm, gamma, uniform, dirichlet
import gaussian_mixture_constr
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
    
  def remove_outliers(self, qtop=0.0025, qbot=0.9975):
    cn2=np.sort(self.data)
    nsamp=cn2.size
    firstel=int(np.floor(nsamp*0.0025))
    lastel=int(np.floor(nsamp*0.9975))
    q_002=cn2[firstel]
    q_998=cn2[lastel]
    if np.abs(q_002-cn2[0])<0.1:
      firstel=0
    if np.abs(q_998-cn2[nsamp-1])<0.1:
      lastel=nsamp-1
    cn2=cn2[firstel:(lastel+1)]
    #print(str(self.comp_id)+"\t"+str(firstel)+"\t"+str(lastel))
    cn2=cn2+norm.rvs(loc=0, scale=0.005, size=cn2.size)
    self.nsamp_rm_outliers=cn2.size
    self.procdata=cn2.reshape(-1, 1)

  def fit_all_models(self, ncarriers):
    fits=[]
    bic_col=6
    mm=1+int(2*np.max(self.procdata))
    if mm<self.nocl_max:
      self.nocl_max=mm
    #nocl_range=np.arange(1, self.nocl_max+1, dtype='int32')
    res=self.fit_one_model(1)
    fits.append(res)
    bic_min=res[bic_col]
    cl_min=1
    nocl=2
    bic=res[bic_col]
    print(str(nocl)+"\t"+str(bic)+"\t"+str(bic_min))
    while (nocl<self.nocl_max) and (bic<=bic_min):
      res=self.fit_one_model(nocl)
      fits.append(res)
      bic=res[bic_col]
      print(str(nocl)+"\t"+str(bic)+"\t"+str(bic_min))
      if bic<bic_min:
        bic_min=bic
        nocl=nocl+1
    cl_min=nocl-1
    print(str(cl_min))
    fits=np.vstack(fits)
    fits=np.hstack([fits, np.empty([fits.shape[0], 3], dtype='float64')])
    fits[:,12]=ncarriers
    fits[:,13]=bic_min
    fits[:,14]=cl_min
    return fits

  def fit_one_model(self, nocl):
    gmm=gaussian_mixture_constr.GaussianMixture(n_components=nocl, n_init=10) 
    gmm.fit(self.procdata)
    mns=gmm.means_
    bic=gmm.bic(self.procdata)
    aic=gmm.aic(self.procdata)
    score=gmm.score(self.procdata)
    kk, mm = gmm.get_kk_mm()
    parity=0
    for ii in range(nocl):
      if ii%2==0:
        parity=parity+gmm.weights_[ii]
    ret=np.array([self.comp_id, self.clus_id, self.clus_dist_id, self.start, self.stop, nocl, bic, aic, score, kk, mm, parity])
    return ret

    

