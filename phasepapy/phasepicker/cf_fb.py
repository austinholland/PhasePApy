import numpy as np
from .util import *
from obspy.signal.filter import bandpass
from obspy.signal.invsim import cosine_taper

class Ratio():
  
  def __init__(self, trace, t_long, freqmin, cnr, perc_taper, mode = 'rms'):
    self.tr = trace
    self.t_long = t_long
    self.freqmin = freqmin
    self.statistics_mode = mode 
    self.cnr = cnr   
    self.perc_taper = perc_taper
    
  def _N_bands(self):
    """ Determine number of band n_bands in term of sampling rate.
    """
    df = self.tr.stats.sampling_rate
    Nyquist = df / 2.0
    n_bands = int(np.log2(Nyquist / 1.5 / self.freqmin)) + 1
    return n_bands

  def filter(self):
    """ Filter data for each band.
    """
    n_bands = self._N_bands()
    LEN = self.tr.stats.npts
    df = self.tr.stats.sampling_rate

    # create zeros 2D array for BF
    BF = np.zeros(shape = (n_bands,LEN))

    for j in range(n_bands):
      octave_high = (self.freqmin + self.freqmin * 2.0) / 2.0 * (2**j)
      octave_low = octave_high / 2.0
      BF[j] = bandpass(self.tr.data, octave_low, octave_high, df, corners = self.cnr, zerophase = False)
      BF[j] = cosine_taper(LEN, self.perc_taper) * BF[j]

    return BF
    
  def _statistics(self):
    """ Calculate statistics for each band.
    """
    n_bands = self._N_bands()
    LEN = self.tr.stats.npts
    dt = self.tr.stats.delta
    npts_t_long = int(self.t_long / dt)
    
    # BF: band filtered data
    BF = self.filter()

    # E: the instantaneous energy
    E = np.power(BF, 2)

    # create zeros 2D array for rmsE, aveE and sigmaE
    rmsE = np.zeros(shape=(n_bands,LEN))
    aveE = np.zeros(shape=(n_bands,LEN))
    sigmaE = np.zeros(shape=(n_bands,LEN))
  
    # range start from 1,not 0, because sigmaE[i=0] calculation encontoured invalid value
    for i in np.arange(1,npts_t_long): 
      rmsE[:,i]=rms(E[:,:i],axis=1)
      aveE[:,i]=np.mean(E[:,:i],axis=1)
      sigmaE[:,i]=np.std(E[:,:i],axis=1)

    # call rolling_window to obtain the rmsE or sigmaE
    tmp=rolling_window(E[:,0:LEN-1],npts_t_long)
    if self.statistics_mode=='rms':
      rmsE[:,npts_t_long:LEN]=rms(tmp,axis=-1)
      rmsE=np.clip(rmsE,1.e-19,1.e19) # clip the sigmaE to avoid zeros values since FC encounters invalid value if denominator sigmaE is zero
    if self.statistics_mode=='std':
      sigmaE[:,npts_t_long:LEN]=np.std(tmp,axis=-1)
      sigmaE=np.clip(sigmaE,1.e-19,1.e19) # clip the sigmaE to avoid zeros values since FC encounters invalid value if denominator sigmaE is zero
      aveE[:,npts_t_long:LEN]=np.mean(tmp,axis=-1)
    
    # calculate statistics
    if self.statistics_mode=='rms':
      FC=np.abs(E)/(rmsE)
    if self.statistics_mode=='std':
      FC=np.abs(E-aveE)/(sigmaE)
    # reassign FC values for the very beginning couple samples to avoid unreasonable large FC from poor sigmaE
    S=self.t_long;L=int(round(S/dt,0))   # S = 0.1
    for k in range(L):
      FC[:,k]=0
    return FC