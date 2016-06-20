""" example_3dassociator_OK.py
This example shows an example using the 1D associator.  The traveltime lookup table for 
this example is included in the data directory.  Third party programs are required to build
the traveltime lookup table.

The associator also uses sqlalchemy to store phase picks and associated events in a database.
This example uses sqlite which should be included in python. 

"""
import sys
sys.path.append("../")
from phasepapy.phasepicker import fbpicker
from phasepapy.associator import tables3D, assoc3D, plot3D
from obspy.core import *
from sqlalchemy.orm import *
from sqlalchemy import create_engine
from datetime import datetime
import glob
import os
# Get logging information
import logging
rootlog=logging.getLogger()
rootlog.setLevel(logging.INFO)
ch=logging.StreamHandler(sys.stderr)
rootlog.addHandler(ch)

# If the associator database exists delete it first, start fresh for this example
if os.path.exists('3dassociator_ok.db'):
  os.remove('3dassociator_ok.db')

# Our SQLite databases are:
db_assoc='sqlite:///3dassociator_ok.db'
db_tt='sqlite:///data_20130616153750/tt_stations_3D.db' # Traveltime database


# Connect to our databases
engine_assoc=create_engine(db_assoc, echo=False)
# Create the tables required to run the 1D associator
tables3D.Base.metadata.create_all(engine_assoc)
Session=sessionmaker(bind=engine_assoc)
session=Session()

# Find all waveform data in the data directory
file_list=glob.glob('data_20130616153750/*.msd')

# Define our picker instance
picker = fbpicker.FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 6, \
    t_up = 0.78, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
 
# Pick the waveforms     
for wf_file in file_list:
  st=read(wf_file)
  #st.merge()  # merge will cause issues if there is a data gap
  for tr in st:
    tr.detrend('linear')
    scnl,picks,polarity,snr,uncert=picker.picks(tr)
    t_create=datetime.utcnow() # Record the time we made the picks
    # Add each pick to the database
    for i in range(len(picks)):
      new_pick=tables3D.Pick(scnl,picks[i].datetime,polarity[i],snr[i],uncert[i],t_create)
      session.add(new_pick) # Add pick i to the database
    session.commit() # Commit the pick to the database
    
# associate picks with phase types
chen3D = assoc3D.LocalAssociator(db_assoc, db_tt, max_km = 350, aggregation = 1, aggr_norm = 'L2', assoc_ot_uncert = 3, nsta_declare = 3, nt = 31, np = 41, nr = 5)
chen3D.id_candidate_events()
chen3D.associate_candidates()

# Plot
plt = plot3D.Plot(db_assoc, db_tt)
plt.cluster_plot(assoc_ot_uncert = 3)
plt.event_plot(1)
plt.section_plot(1, files)