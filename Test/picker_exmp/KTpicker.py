#!/usr/bin/env python

from PhasePApy.PhasePicker import ktpicker
from obspy.core import read

#=======================================================================
# KTpicker example
file = open('../20130616153750/OKCFA.HHZ.OK...20130616153750.msd')
st = read(file); tr = st[0]; tr.detrend('linear')

chenPicker = ktpicker.KTPicker(t_win = 1, t_ma = 10, nsigma = 6, t_up = 0.78, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
scnl, picks, polarity, snr, uncert = chenPicker.picks(tr)
print 'scnl:', scnl
print 'picks:', picks 
print 'polarity:', polarity
print 'signal to noise ratio:', snr
print 'uncertainty:', uncert

# plot
summary = ktpicker.KTSummary(chenPicker, tr)
summary.plot_summary()
summary.plot_picks()
