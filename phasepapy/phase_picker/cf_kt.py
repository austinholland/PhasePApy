import numpy as np
from scipy.stats import kurtosis
from .util import *

class Kurtosis():
  
  def __init__(self, trace, t_win):
    self.tr = trace
    self.t_win = t_win
    self.npts = self.tr.stats.npts
    self.sampling_rate = self.tr.stats.sampling_rate
    self.delta = 1.0/self.tr.stats.sampling_rate
  
  def _statistics(self):
    data = self.tr.data
    t = np.arange(0, self.delta * self.npts, self.delta)
    m = len(data)
    Nsta = int(self.t_win * self.sampling_rate)

    # compute the short time average (STA)
    kt = np.zeros(m, dtype='float64')
    pad_kt = np.zeros(Nsta)
    # Tricky: Construct a big window of length len(a)-nsta. Now move this
    # window nsta points, i.e. the window "sees" every point in a at least
    # once.
    # Changed xrange to range as it is compatible in both python 2 & 3
    for i in range(m):  # window size to smooth over
        kt[i] = abs(kurtosis(data[i-Nsta:i]))

    kt[0:Nsta] = 0
    
    return kt