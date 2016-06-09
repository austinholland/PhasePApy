""" example_ktpicker.py

Usage: python example_ktpicker.py
"""
import sys
sys.path.append("../")
from phasepapy.phasepicker import ktpicker
from obspy.core import *
import obspy.clients.iris as iris

#=======================================================================
# KTpicker example
# Load data into an Obspy Stream
# Obspy version dependent check the documentation for your version
wfstart=UTCDateTime(2016,6,9,11,11,28)
iris_client=iris.Client()

st=iris_client.timeseries("OK","CROK","--","HHZ",wfstart,wfstart+10*60) # Get ten minutes of data
st.merge() # Ensure that traces aren't split
tr=st[0]
tr.detrend('linear') # Perform a linear detrend on the data

chenPicker = ktpicker.KTPicker(t_win = 1, t_ma = 10, nsigma = 6, t_up = 0.78, nr_len = 2, 
  nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
scnl, picks, polarity, snr, uncert = chenPicker.picks(tr)
print('scnl:', scnl)
print('picks:', picks)
print('polarity:', polarity)
print('signal to noise ratio:', snr)
print('uncertainty:', uncert)

# plot
summary = ktpicker.KTSummary(chenPicker, tr)
summary.plot_summary()
summary.plot_picks()
