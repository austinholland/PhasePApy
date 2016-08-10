""" full_example1D.py

This example provides an example of making a travel-time lookup table and picking data from
one single source using built-in obspy functionality.  The example is actually written
such that some of the methods could be moved elsewhere or called by other python scripts.

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
import obspy.taup as taup
import numpy as np
import re
import os


def build_tt_tables(minlat=None,maxlat=None,minlon=None,maxlon=None,channel_codes=['EH','BH','HH'],db=None,maxdist=500.,source_depth=5.):
  """ channel_codes select channels that start with those codes  
  maximum distance is in km
  source depth is generally set to the average earthquake depth for the region you are working
  for more granularity use the 3D associator
  """
  # Create a connection to an sqlalchemy database
  tt_engine=create_engine(db,echo=False)
  tt_stations_1D.Base.metadata.create_all(tt_engine)
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
        print(ch)
        print(dir(ch))
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
  return True
      
  





if __name__=="__main__":
  # Run this simple example
  # If the associator database exists delete it first, start fresh for this example
  if os.path.exists('1dassociator_ex.db'):
    os.remove('1dassociator_ex.db')

  # Our SQLite databases are:
  db_assoc='sqlite:///1dassociator_ex.db'
  db_tt='sqlite:///tt_ex_1D.db' # Traveltime database


  # Connect to our databases
  engine_assoc=create_engine(db_assoc, echo=False)
  # Create the tables required to run the 1D associator
  tables1D.Base.metadata.create_all(engine_assoc)
  Session=sessionmaker(bind=engine_assoc)
  session=Session()

  build_tt_tables(minlat=34.0,maxlat=37.2,minlon=-101.0,maxlon=-95.0,channel_codes=['EH','BH','HH'],db=db_tt)