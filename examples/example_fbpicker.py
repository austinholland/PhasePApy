""" example_fbpicker.py
Example using the FBpicker to pick one channel of data examples for a 3-component station.
Data is downloaded from IRIS DMC using obspy check your version of obspy for compatability.

Usage:  python example_fbpicker.py
"""
import sys
sys.path.append("../")
from phasepapy.phasepicker import fbpicker
from obspy.core import *
import obspy.clients.iris as iris

# Load data into an Obspy Stream
# Obspy version dependent check the documentation for your version
wfstart=UTCDateTime(2016,6,9,11,11,28)
iris_client=iris.Client()

# Example on the vertical component
st=iris_client.timeseries("OK","CROK","--","HHZ",wfstart,wfstart+10*60) # Get ten minutes of data
st.merge() # Ensure that traces aren't split
tr=st[0]
tr.detrend('linear') # Perform a linear detrend on the data

# Create a picker object with pick parameters
picker = fbpicker.FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 8, \
  t_up = 0.4, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
scnl, picks, polarity, snr, uncert = picker.picks(tr) 
print('scnl:', scnl)
print('picks:', picks)
print('polarity:', polarity)
print('signal to noise ratio:', snr)
print('uncertainty:', uncert)

# plot
summary = fbpicker.FBSummary(picker, tr)
summary.plot_bandfilter()
summary.plot_statistics()
summary.plot_summary()

# Example for horizontal components

st=iris_client.timeseries("OK","CROK","--","HHN",wfstart,wfstart+10*60) # Get ten minutes of data
st.merge() # Ensure that traces aren't split
tr=st[0]
tr.detrend('linear') # Perform a linear detrend on the data

# Create a picker object with pick parameters
picker = fbpicker.FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 8, \
  t_up = 0.4, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
scnl, picks, polarity, snr, uncert = picker.picks(tr) 
print('scnl:', scnl)
print('picks:', picks)
print('polarity:', polarity)
print('signal to noise ratio:', snr)
print('uncertainty:', uncert)

# plot
summary = fbpicker.FBSummary(picker, tr)
summary.plot_bandfilter()
summary.plot_statistics()
summary.plot_summary()

st=iris_client.timeseries("OK","CROK","--","HHE",wfstart,wfstart+10*60) # Get ten minutes of data
st.merge() # Ensure that traces aren't split
tr=st[0]
tr.detrend('linear') # Perform a linear detrend on the data

# Create a picker object with pick parameters
picker = fbpicker.FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 8, \
  t_up = 0.4, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
scnl, picks, polarity, snr, uncert = picker.picks(tr) 
print('scnl:', scnl)
print('picks:', picks)
print('polarity:', polarity)
print('signal to noise ratio:', snr)
print('uncertainty:', uncert)

# plot
summary = fbpicker.FBSummary(picker, tr)
summary.plot_bandfilter()
summary.plot_statistics()
summary.plot_summary()