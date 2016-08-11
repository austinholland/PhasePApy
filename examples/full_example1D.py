""" full_example1D.py

This example provides an example of making a travel-time lookup table and picking data from
one single source using built-in obspy functionality.  The example is actually written
such that some of the methods could be moved elsewhere or called by other python scripts.

As a warning if you go through and pick and associate the data this script takes a while, 
but it does demonstrate how to build a 1D velocity model look up table and then pick the
data.

Austin Holland
"""
import sys
sys.path.append("../")
from phasepapy.phasepicker import fbpicker
from phasepapy.associator import tables1D, assoc1D, plot1D
from phasepapy.associator import tt_stations_1D
from sqlalchemy.orm import *
from sqlalchemy import create_engine
from obspy.core import *
import obspy.geodetics as geodetics
import obspy.clients.fdsn as fdsn
import obspy.clients.iris as iris
import obspy.taup as taup
import numpy as np
from datetime import datetime
import re
import os

# Get logging information
import logging
rootlog=logging.getLogger()
rootlog.setLevel(logging.INFO)
ch=logging.StreamHandler(sys.stderr)
rootlog.addHandler(ch)

def build_tt_tables(minlat=None,maxlat=None,minlon=None,maxlon=None,channel_codes=['EH','BH','HH'],db=None,maxdist=500.,source_depth=5.):
  """ channel_codes select channels that start with those codes  
  maximum distance is in km
  source depth is generally set to the average earthquake depth for the region you are working
  for more granularity use the 3D associator
  """
  # Create a connection to an sqlalchemy database
  tt_engine=create_engine(db,echo=False)
  tt_stations_1D.BaseTT1D.metadata.create_all(tt_engine)
  TTSession=sessionmaker(bind=tt_engine)
  tt_session=TTSession()
  # Create a cliet to IRIS FDSN
  fdsnclient=fdsn.Client("IRIS")
  # Create an obspy inventory of stations
  #http://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_stations.html#obspy.clients.fdsn.client.Client.get_stations
  inv=fdsnclient.get_stations(minlatitude=minlat,maxlatitude=maxlat,minlongitude=minlon,maxlongitude=maxlon,level='channel')
  # Plot our results just for fun
  inv.plot(projection='ortho',color_per_network='True')
  # Now save these station into the 1D travel-time table database
  # The associator could be modified to interact with Obspy Inventory objects
  for net in inv:
    network=net.code
    for sta in net:
      loccodes=[]
      for ch in sta:
#         print(ch)
#         print(dir(ch))
        for cc in channel_codes:
          if re.match(cc,ch.code):
            if not ch.location_code in loccodes:
              loccodes.append(ch.location_code)
      for loc in loccodes:  
        station=tt_stations_1D.Station1D(sta.code,network,loc,sta.latitude,sta.longitude,sta.elevation)
        # Save the station locations in the database
        tt_session.add(station)
      tt_session.commit()

  # Now we have to build our traveltime lookup tables
  # We will use IASP91 here but obspy.taup does let you build your own model
  velmod=taup.TauPyModel(model='iasp91')      
  # Define our distances we want to use in our lookup table
  delta_distance=1. # km for spacing tt calculations  
  # Probably better to use a progressive type scheme instead of linear, but this is an example
  distance_km=np.arange(0,maxdist+delta_distance,delta_distance)
  for d_km in distance_km:
    d_deg=geodetics.kilometer2degrees(d_km)
    ptimes=[]
    stimes=[]
    p_arrivals=velmod.get_travel_times(source_depth_in_km=source_depth,
      distance_in_degree=d_deg,phase_list=['P','p'])
    for p in p_arrivals:
      ptimes.append(p.time)
    s_arrivals=velmod.get_travel_times(source_depth_in_km=source_depth,
      distance_in_degree=d_deg,phase_list=['S','s'])
    for s in s_arrivals:
      stimes.append(s.time)
    tt_entry=tt_stations_1D.TTtable1D(d_km,d_deg,np.min(ptimes),np.min(stimes),np.min(stimes)-np.min(ptimes))
    tt_session.add(tt_entry)
    tt_session.commit() # Probably faster to do the commit outside of loop but oh well
  tt_session.close()
  return inv
      
def irisws_pick(dbsession=None,picker=None,inventory=None,starttime=None,endtime=None,t_chunk=1200,channel_codes=['EH','BH','HH'],overlap=30.):
  """  Each time chunk is overlapped some and so is technically longer than each t_chunk by that amount"""  
  irisclient=iris.Client()
  inv=inventory.select(starttime=starttime,endtime=endtime)
  for net in inv:
    network=net.code
    for sta in net:
      for ch in sta:
        for cc in channel_codes:
          if re.match(cc,ch.code):
            tlast=starttime
            while tlast<endtime:  #Loop through all the time we want to pick
              if tlast+t_chunk>endtime:
                end=endtime
              else:
                end=tlast+t_chunk
                #We decimate the data as this is an example this will keep memory usage down
                try:
                  st=irisclient.timeseries(network,sta.code,ch.location_code,ch.code,tlast-30,end)
                  logging.info(str(st))
                except:
                  logging.info("%s.%s.%s.%s not available" % (network,sta.code,ch.location_code,ch.code))
                if len(st)>0:  #Make sure we got data
                  for tr in st:  # For each data segment run the picker
                    tr.detrend('linear')
                    scnl,picks,polarity,snr,uncert=picker.picks(tr)
                    t_create=datetime.utcnow() # Record the time we made the picks
                    # Add each pick to the database
                    for i in range(len(picks)):
                      new_pick=tables1D.Pick(scnl,picks[i].datetime,polarity[i],snr[i],uncert[i],t_create)
                      dbsession.add(new_pick) # Add pick i to the database
                    dbsession.commit() # Commit the pick to the database

              tlast=end
            




if __name__=="__main__":
  # Run this simple example
  # If the associator database exists delete it first, start fresh for this example
  if os.path.exists('1dassociator_ex.db'):
    os.remove('1dassociator_ex.db')
    os.remove('tt_ex_1D.db')

  # Our SQLite databases are:
  db_assoc='sqlite:///1dassociator_ex.db'
  db_tt='sqlite:///tt_ex_1D.db' # Traveltime database


  # Connect to our databases
  engine_assoc=create_engine(db_assoc, echo=False)
  # Create the tables required to run the 1D associator
  tables1D.Base.metadata.create_all(engine_assoc)
  Session=sessionmaker(bind=engine_assoc)
  session=Session()

  inventory=build_tt_tables(minlat=34.0,maxlat=37.2,minlon=-101.0,maxlon=-95.0,channel_codes=['EH','BH','HH'],db=db_tt)
  
  # Define our picker instance
  picker = fbpicker.FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 7, \
    t_up = 0.7, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
  start=UTCDateTime('2016-08-09T21:40:00.0')
  end=UTCDateTime('2016-08-09T22:20:00.0')  
  
  irisws_pick(dbsession=session,
    picker=picker,
    inventory=inventory,
    starttime=start,
    endtime=end,
    t_chunk=1200,
    channel_codes=['EH','BH','HH'],
    overlap=30.)
  
  # Define the associator
  assocOK=assoc1D.LocalAssociator(db_assoc, db_tt,
    max_km = 350, 
    aggregation = 1, 
    aggr_norm = 'L2', 
    cutoff_outlier = 10, 
    assoc_ot_uncert = 7, 
    nsta_declare = 4, 
    loc_uncert_thresh = 0.2)
  
  # Identify candidate events (Pick Aggregation)
  assocOK.id_candidate_events()

  # Associate events
  assocOK.associate_candidates()

  # Add singles stations to events
  assocOK.single_phase()

  # Plot example event
  plt=plot1D.Plot(db_assoc,db_tt)
  plt.cluster_plot(assoc_ot_uncert=3)
  plt.event_plot(1,west = -104.5, east= -94, south = 33.5, north = 37.5, deltalon = 1.0, deltalat = 1.0)
  plt.section_plot(1,file_list,seconds_ahead = 5, record_length = 100, channel = 'Z')
  