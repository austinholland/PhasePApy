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
from phasepapy.associator.tt_stations_1D import TTtable1D, Station1D
from sqlalchemy.orm import *
from sqlalchemy import create_engine
from obspy.core import *
import obspy.clients.fdsn as fdsn
import obspy.taup as taup
import re


def build_tt_tables(minlat=None,maxlat=None,minlon=None,maxlon=None,channel_codes=['EH','BH','HH'],db=None,maxdist=500.):
  """ channel_codes select channels that start with those codes  
  maximum distance is in km
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
        print(ch)
        print(dir(ch))
        for cc in channel_codes:
          if re.match(cc,ch.code):
            if not ch.location_code in loccodes:
              loccodes.append(ch.location_code)
      for loc in loccodes:  
        station=Station1D(sta.code,network,loc,sta.latitude,sta.longitude,sta.elevation)
        # Save the station locations in the database
        tt_session.add(station)
      tt_session.commit()

  # Now we have to build our traveltime lookup tables
  velmod=taup.TauPyModel(model='iasp91')      
  # Define our distances we want to use in our lookup table
  delta_d=1. # km for spacing tt calculations  
  # Probably better to use a progressive type scheme instead of linear, but this is an example
  d_km=np.arange(0,maxdist+delta_d,delta_d)
  
      
  





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