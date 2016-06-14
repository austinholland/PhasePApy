#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib
import copy
from obspy.core import *
from .scnl import * 
from .cf_fb import *


class FBPicker():
  """ 
  FBpicker is designed to make picks based on the energy ratio of filtered data.

  This algorithm is scripted by Chen Chen and Guided by Austin A. Holland, June 2015.

  This is a pure python implementation of FBpicker it might be more efficient to wrap the 
  original c-code with swig, but it appears to be more complex than necessary. It uses 
  obspy and numpy to implement the picker routines. 

  FBpicker could be used separately just for picking picks. Also the Associator can 
  follow to associate the picks with phase types.

  """
  
  def __init__(self, t_long = 5, freqmin = 1, corner = 1, perc_taper = 0.1, mode = 'rms', t_ma = 20, nsigma = 6,\
      t_up = 0.78, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3):
    
    """ 
    Parameter description:
      t_long       : the time in seconds of moving window to calculate CFn of each bandpass filtered data
      freq_min     : the center frequency of first octave filtering band
      cnr          : corner order of bandpass filtering
      t_ma         : the time in seconds of the moving average window for dynamic threshold
      n_sigma      : controls the level of threshold to trigger potential picks
      t_up         : the time in seconds not allowed consecutive pick in this duration
      mode         : two options: standard deviation (std) or root mean square (rms)
      nr_len       : noise ratio filter window length before and after potential picks used to calculate standard deviation 
      nr_coeff     : control threshold level to determine if remove the pick by comparing std or rms on both sides of each potential pick  
      pol_len      : window length in samples to calculate the standard deviation of waveform before the picks
      pol_coeff    : determine if declare first motion as 'Compression' or 'Dilation' by comparing the first local extreme value after pick and standard deviation in previous window
      uncert_len   : window length in time to calculate the rms of the CF before the picks, we make it as long as t_ma
      uncert_coeff : control the floating level based on the noise of CF
    """

    self.t_long = t_long
    self.freqmin = freqmin
    self.cnr = corner
    self.perc_taper = perc_taper
    self.statistics_mode = mode
    self.t_ma = t_ma
    self.nsigma = nsigma
    self.t_up = t_up
    self.nr_len = nr_len  
    self.nr_coeff = nr_coeff
    self.pol_len = pol_len 
    self.pol_coeff = pol_coeff
    self.uncert_len = self.t_ma
    self.uncert_coeff = uncert_coeff
    

  def picks(self,tr):
    """ 
    Make picks, polarity, snr, and uncertainty.
    """
    
    #tr = trace.detrend('linear')
#       now=time.time()
    summary=FBSummary(self,tr)
#       print "It took %f s to summary" %(time.time()-now)

    # threshold
#       now=time.time()
    threshold=summary.threshold
#       print "It took %f s to threshold" %(time.time()-now)

    # picks
#       now=time.time()
    scnl,picks,trigger,snr=summary.pick_ident()
#       print "It took %f s to ident" %(time.time()-now)

    # uncertainty
#       now=time.time()
    uncertainty=summary.uncertainty()
#       print "It took %f s to uncertainty" %(time.time()-now)
  
    # polarity
#       now=time.time()
    polarity=summary.polarity()
#       print "It took %f s to polarity" %(time.time()-now)
    return scnl,picks,polarity,snr,uncertainty  



class FBSummary():
  """ 
  The class summarizes all bands CFn, threshold level, cleans the false picks, determines uncertainty, polarity 
  and plots band filtered data, statistics for each band, and CF.
  """
  
  def __init__(self,picker, tr):
    self.picker = picker
    self.tr = tr
    self.stats = self.tr.stats
    self.t_long = self.picker.t_long
    self.cnr = self.picker.cnr
    self.perc_taper = self.picker.perc_taper
    self.cf = Ratio(self.tr, self.picker.t_long, self.picker.freqmin, self.cnr, self.perc_taper, self.picker.statistics_mode)
    self.FC = self.cf._statistics()
    self.summary = np.amax(self.FC,axis=0)
    self.thres = self.threshold() 
    self.uncert = self.uncertainty()
    self.pol = self.polarity()  
    
  
  def threshold(self):
    """ 
    Control the threshold level with nsigma.
    """
    dt = self.stats.delta
    npts_Tma = int(round(self.picker.t_ma / dt,0))
    LEN = self.stats.npts
    #print "npts_Tma: ",npts_Tma
    threshold = np.zeros(LEN)
    #threshold[0:npts_Tma] = 1
    threshold[npts_Tma:LEN] = rms(rolling_window(self.summary[0:LEN - 1],npts_Tma), -1) * self.picker.nsigma

    return threshold
  
  def pick_ident(self):
    """ 
    Clean false picks and Make picks.
    """
    
    scnl = SCNL([self.stats.station,self.stats.channel,self.stats.network,self.stats.location])
    dt = self.stats.delta
    npts_Tma = int(round(self.picker.t_ma/dt,0))
    LEN = self.stats.npts
    
    # trigger the earthquakes
    trigger_ptnl_index = np.where(self.summary[npts_Tma:LEN] > self.thres[npts_Tma:LEN])
    trigger_ptnl_index = trigger_ptnl_index + np.array(npts_Tma)
    t = np.arange(0, self.stats.npts/self.stats.sampling_rate, dt) 
    trigger_ptnl = t[trigger_ptnl_index][0]
    
    # clean close picks
    window_t_up = int(round(self.picker.t_up/dt,0))
    trigger_remove1_index = []	
    for i in range(0,len(trigger_ptnl) - 1): # second from last
      if (trigger_ptnl[i+1] - trigger_ptnl[i]) <= window_t_up * dt: # avoid consecutive picking
        trigger_remove1_index.append(i+1)
    # delete close picks
    trigger_ptnl=np.delete(trigger_ptnl, trigger_remove1_index)
    
    # clean_filtering
    trigger_remove2_index = []
    N = self.picker.nr_coeff
    filter_length = self.picker.nr_len    
    for i in range(len(trigger_ptnl)):  
      # determine filter_length for each pick:
      r, R = self.winlen(i, trigger_ptnl, filter_length, t, dt)
      M = min(r, R)
      if N * np.std(self.tr.data[int(round(trigger_ptnl[i] / dt,0)) - M:int(round(trigger_ptnl[i] / dt,0))]) >= np.std(self.tr[int(round(trigger_ptnl[i] / dt,0)):int(round(trigger_ptnl[i] / dt,0)) + M]): 
        trigger_remove2_index.append(i)  
    #delete false picks	
    trigger_ptnl = np.delete(trigger_ptnl, trigger_remove2_index)
    
    # assign potential picks to trigger
    trigger = trigger_ptnl
    
    # really need to be careful copy list or array, since if just copy like A=B, when any element in B changes, A will change as well
    # roll backward for picking
    picks = []
    for i in range(len(trigger)):
      index = int(round(trigger[i] / dt,0))
      while True:
        if self.summary[index] > self.summary[index - 1]:
          index -= 1
        else:
          break
      picks.append(UTCDateTime(self.tr.stats.starttime + round(t[index], 3)))

    # really need to be careful copy list or array, since if just copy like A=B, when any element in B changes, A will change as well
    # roll forward for maximum signal values
    maxes = copy.deepcopy(trigger)
    for i in range(len(trigger)):
       index = int(round(trigger[i] / dt,0))
       while True:
         if self.summary[index] < self.summary[index + 1]:
           index += 1
         else:
           break
       maxes[i] = round(self.summary[index],3)
    
    # really need to be careful copy list or array, since if just copy like A=B, when any element in B changes, A will change as well   
    # Signal noise ration: SNR
    SNR=copy.deepcopy(trigger)

    for i in range(len(picks)):
      index = int(round(trigger[i] / dt,0))
      noise = rms(self.summary[index-npts_Tma:index])
      SNR[i] = round(maxes[i] / noise,1)
    
    return scnl,picks,trigger,SNR
  
  def uncertainty(self):
    """ 
    Uncertainty is determined based on the noise level of CF.
    """
    
    scnl,picks,trigger,SNR = self.pick_ident()
    dt = self.stats.delta
    npts_Tma = int(round(self.picker.t_ma / dt, 0))
    t = np.arange(0, self.stats.npts / self.stats.sampling_rate, dt) 
    pick_uncert = copy.deepcopy(trigger)

    for i in range(len(trigger)):
      r, R=self.winlen(i, trigger, npts_Tma, t, dt)
      index0 = int(round((picks[i] - self.tr.stats.starttime) / dt, 0))
      index = int(round(trigger[i] / dt, 0))
      uncert_level = self.picker.uncert_coeff * rms(self.summary[index0 - npts_Tma:index0])
      while True:
        if self.summary[index] > uncert_level and self.summary[index] > self.summary[index - 1]:
          index -= 1
        else:
          break
      
      # t[index] is the last sample lower than uncert_level, t[index+1] is the first sample higher than uncert_level
      pick_uncert[i] = round(t[index+1] - (picks[i] - self.tr.stats.starttime), 3)    
         
    return pick_uncert
  
  def polarity(self):
    """ 
    Determine polarity for declared picks.
    """
    dt=self.stats.delta
    t = np.arange(0, self.stats.npts/self.stats.sampling_rate, dt) 
    pol=[]
    scnl,picks,trigger,snr=self.pick_ident()
    for i in range(len(picks)):
      index0=int(round((picks[i]-self.tr.stats.starttime)/dt,0))
      index=index0
      
      # roll forward index+=1
      while True:
        if index>=self.stats.npts-1-2:
          break
        elif (self.tr[index+1]-self.tr[index])*(self.tr[index+2]-self.tr[index+1])>0:
          index+=1
        else:
          break
      
      # notice index+1, rolling stop one point before extreme, compare with std to avoid very small 
      if self.tr[index+1] - self.tr[index0] > 0 and abs(self.tr[index+1] - self.tr[index0]) > self.picker.pol_coeff * np.std(self.tr[index0 - self.picker.pol_len: index0]):
        polarity='C'
      elif self.tr[index+1] - self.tr[index0] < 0 and abs(self.tr[index+1] - self.tr[index0]) > self.picker.pol_coeff * np.std(self.tr[index0 - self.picker.pol_len: index0]):
        polarity='D'
      else: 
        polarity=''
      pol.append(polarity)     
    return pol
  
  def winlen(self,index,trigger_ptnl,filter_length,t,dt):
    """ 
    Determine the filter window length. If the time difference between two picks is less 
    than window length, use the picks interval as window.
    """
    i=index
    if len(trigger_ptnl)==1:
      # print 'A'
      if trigger_ptnl[i]<=filter_length:
        r=int(round(trigger_ptnl[i]/dt,0))
      if trigger_ptnl[i]>filter_length:
        r=int(round(filter_length/dt,0))
      if trigger_ptnl[i]+filter_length>=t[-1]: 
        R=int(round((t[-1]-trigger_ptnl[i])/dt,0))
      if trigger_ptnl[i]+filter_length<t[-1]:
        R=int(round(filter_length/dt,0))
    elif len(trigger_ptnl)>1:
      # print 'B'
      if i==0:
        if trigger_ptnl[i]<=filter_length: 
          r=int(round(trigger_ptnl[i]/dt,0)); # print 'a'
        if  trigger_ptnl[i]>filter_length:
          r=int(round(filter_length/dt,0));  # print 'b'
        if (trigger_ptnl[i+1]-trigger_ptnl[i])<=filter_length:
          R=int(round((trigger_ptnl[i+1]-trigger_ptnl[i])/dt,0));  # print 'c'
        if (trigger_ptnl[i+1]-trigger_ptnl[i])>filter_length:
          R=int(round(filter_length/dt,0));  # print 'd'
      elif i>0 and i<len(trigger_ptnl)-1:
        if trigger_ptnl[i]-trigger_ptnl[i-1]<=filter_length:
          r=int(round((trigger_ptnl[i]-trigger_ptnl[i-1])/dt,0));  # print 'e'
        if trigger_ptnl[i]-trigger_ptnl[i-1]>filter_length:
          r=int(round(filter_length/dt,0));  # print 'f'
        if trigger_ptnl[i+1]-trigger_ptnl[i]<=filter_length:
          R=int(round((trigger_ptnl[i+1]-trigger_ptnl[i])/dt,0));  # print 'g'
        if trigger_ptnl[i+1]-trigger_ptnl[i]>filter_length:
          R=int(round(filter_length/dt,0));  # print 'h'
      elif i==len(trigger_ptnl)-1:
        if trigger_ptnl[i]-trigger_ptnl[i-1]<=filter_length:
          r=int(round((trigger_ptnl[i]-trigger_ptnl[i-1])/dt,0));  # print 'i'
        if trigger_ptnl[i]-trigger_ptnl[i-1]>filter_length:
          r=int(round(filter_length/dt,0));  # print 'j'
        if trigger_ptnl[i]+filter_length>t[-1]: 
          R=int(round((t[-1]-trigger_ptnl[i])/dt,0));  # print 'k'
        if trigger_ptnl[i]+filter_length<=t[-1]:
          R=int(round(filter_length/dt,0));  # print 'l'
    return r,R  
    
  def plot_bandfilter(self):
    """ 
    Plot bandpass filtered data BF for each band.
    """
    
    matplotlib.rcParams["axes.labelsize"]="large"
    matplotlib.rcParams["axes.linewidth"]=2.0
    matplotlib.rcParams["xtick.major.size"]=8
    matplotlib.rcParams["ytick.major.size"]=8
    matplotlib.rcParams["ytick.minor.size"]=5
    matplotlib.rcParams["xtick.labelsize"]="large"
    matplotlib.rcParams["ytick.labelsize"]="large"
    
    dt=self.stats.delta
    t = np.arange(0,self.stats.npts/self.stats.sampling_rate,dt)
    n_bands=self.cf._N_bands()
    BF=self.cf.filter()
    fig=plt.figure(figsize=(10,12))
    title='{0} {1} {2}'.format(self.stats.station, self.stats.channel, self.stats.starttime)
    
    # Plot raw data
    ax = plt.subplot(n_bands+1,1,1)
    ax.axes.get_xaxis().set_visible(False)
    ax.plot(t, self.tr, c='gray')
    plt.ylabel('Raw Data')
    plt.title(title)
    
    # Plot filtered data
    st=self.picker.freqmin # for subtitle
    unit='Hz'
    for i in range(n_bands):
        ax = plt.subplot(n_bands+1,1,i+2)
        ax.plot(t, BF[i], c='k')
        subtitle=str(st*2**i)+unit
        ax.text(0.5,0.8,subtitle,horizontalalignment='center',transform=ax.transAxes)
        ax.axes.get_xaxis().set_visible(False)
        plt.ylabel('BF%s'%(i+1))
    ax.axes.get_xaxis().set_visible(True)
    plt.xlabel('Time (s)')
    #fig.savefig('filtering data.pdf')
    #plt.suptitle(title)
    plt.tight_layout()
    plt.show()  
  
  def plot_statistics(self):
    """ 
    Plot CFn for each band.
    """
    
    matplotlib.rcParams["axes.labelsize"]="large"
    matplotlib.rcParams["axes.linewidth"]=2.0
    matplotlib.rcParams["xtick.major.size"]=8
    matplotlib.rcParams["ytick.major.size"]=8
    matplotlib.rcParams["ytick.minor.size"]=5
    matplotlib.rcParams["xtick.labelsize"]="large"
    matplotlib.rcParams["ytick.labelsize"]="large"

    dt = self.stats.delta
    t = np.arange(0, self.stats.npts/self.stats.sampling_rate, dt) 
    scnl,picks,trigger,snr = self.pick_ident()
    fig = plt.figure(figsize=(10,12))
    title='{0} {1} {2}'.format(self.stats.station, self.stats.channel, self.stats.starttime)
    
    # Plot raw data
    n_bands = self.cf._N_bands()
   
    ax = plt.subplot(n_bands+1,1,1)
    ax.plot(t, self.tr, c='gray')
    ax.axes.get_xaxis().set_visible(False)
    # plot picks
    #for i in range(len(picks)):
    #  ax.plot([(picks[i]-self.tr.stats.starttime), (picks[i]-self.tr.stats.starttime)], [min(self.tr),max(self.tr)], 'r--')
    #  ax.text((picks[i]-self.tr.stats.starttime),max(self.tr)-0.3*(max(self.tr)-min(self.tr)),'%s' % (self.pol[i]),color='red')
    plt.ylabel('Raw Data')
    plt.title(title)
    
    # Plot filtered data
    freq=self.picker.freqmin # for subtitle
    unit='Hz'
    for i in range(n_bands):
        ax = plt.subplot(n_bands+1,1,i+2)
        ax.plot(t, self.FC[i], c='k')
        subtitle=str(freq*2**i)+unit
        ax.text(0.5,0.8,subtitle,horizontalalignment='center',transform=ax.transAxes)
        ax.axes.get_xaxis().set_visible(False)
        plt.ylabel('CF%s'%(i+1))
    ax.axes.get_xaxis().set_visible(True)
    plt.xlabel('Time (s)')
    #fig.savefig('CFn_FB.pdf')
    #plt.suptitle(title)
    plt.tight_layout()
    plt.show()
  
  def plot_summary(self):
    """ 
    Plot CF.
    """
    
    matplotlib.rcParams["axes.labelsize"]="large"
    matplotlib.rcParams["axes.linewidth"]=2.0
    matplotlib.rcParams["xtick.major.size"]=8
    matplotlib.rcParams["ytick.major.size"]=8
    matplotlib.rcParams["ytick.minor.size"]=5
    matplotlib.rcParams["xtick.labelsize"]="large"
    matplotlib.rcParams["ytick.labelsize"]="large"
    
    fig=plt.figure(figsize=(12,8))
    dt=self.stats.delta
    t = np.arange(0, self.stats.npts/self.stats.sampling_rate, dt)  
    
    # Plot raw data
    ax = plt.subplot(2,1,1)
    ax.plot(t, self.tr, c='gray')
    plt.ylabel('Raw Data')
    
    # Plot summary
    ax1 = plt.subplot(2,1,2)
    ax1.plot(t,self.summary/np.amax(self.summary),c='k')
    ax1.plot(t,self.thres/np.amax(self.summary),'--',linewidth=2.0,c='k')
    ax1.legend(('Normalized CF','Threshold','Picks'),'upper right', shadow=True, fancybox=True)
    plt.ylabel('Characteristic Function')
    plt.xlabel('Time (s)')
    #plt.xlim([65,69])   # zoom in to see the polarity
    
    # Plot picks
    scnl,picks,trigger,snr=self.pick_ident()
    for i in range(len(picks)):
      ax.plot([(picks[i]-self.tr.stats.starttime),(picks[i]-self.tr.stats.starttime)],[min(self.tr.data),max(self.tr.data)],'k--')
      ax.text((picks[i]-self.tr.stats.starttime),1,'%s' % (self.pol[i]),color='black')
      
    plt.tight_layout()
    plt.show()

