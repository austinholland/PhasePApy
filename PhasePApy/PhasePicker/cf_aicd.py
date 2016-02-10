import numpy as np
import math
from util import *

class AicDeriv():
  
  def __init__(self, trace):
    self.tr = trace  
    
  def _statistics(self):
    npts = self.tr.stats.npts
    data = self.tr.data;
    delta = 1.0/self.tr.stats.sampling_rate
    AIC = np.zeros(npts)
    
    # reverse indexing to remove the nan, np.std(data[:0]) is nan, starting index need to 
    # be npts-2, if data array only has 1 sample, the std is 0, the log10(0) is inf
    for k in range(npts-2,0,-1):

      a = k*np.log10(np.std(data[:k])**2)+(npts-k-1)*np.log10(np.std(data[k:])**2)
      
      #print a,np.log10(np.std(data[k:]))
      if a == -float('inf'):
        a = AIC[k+1]
      AIC[k] = a
    AIC[0] = AIC [1]
    AIC[-1] = AIC[-2]
    AIC = np.array(AIC)
    
    AIC_deriv = []
    for i in range(npts-1):
      b = np.abs(AIC[i+1]-AIC[i])
      AIC_deriv.append(b)

    AIC_deriv.insert(0,AIC_deriv[0])
    AIC_deriv = np.array(AIC_deriv)

    return AIC, AIC_deriv