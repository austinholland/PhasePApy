#!/usr/bin/env python

from PhasePApy.PhasePicker import fbpicker
from obspy.core import read

#=======================================================================
# FBpicker example
file = open('../20130616153750/OKCFA.HHZ.OK...20130616153750.msd')
st = read(file); tr=st[0]; tr.detrend('linear')

chenPicker = fbpicker.FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 6, \
  t_up = 0.4, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
scnl, picks, polarity, snr, uncert = chenPicker.picks(tr) 
print 'scnl:', scnl
print 'picks:', picks 
print 'polarity:', polarity
print 'signal to noise ratio:', snr
print 'uncertainty:', uncert

# plot
summary = fbpicker.FBSummary(chenPicker, tr)
summary.plot_bandfilter()
summary.plot_statistics()
summary.plot_summary()
