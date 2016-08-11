from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base

BaseTT1D=declarative_base()

class TTtable1D(BaseTT1D):
  __tablename__='traveltimes'
  id=Column(Integer, primary_key=True)
  d_km=Column(Float)     # Distance in kilometers
  delta=Column(Float)    # Distance in degrees
  p_tt=Column(Float)     # P travel-time
  s_tt=Column(Float)     # S travel-time
  s_p=Column(Float)      # S-P time
  
  def __init__(self,d_km,delta,p_tt,s_tt,s_p):
    """ TTtable(d_km,delta,p_tt,s_tt,s_p)
    Create an entry in the travel-time table for later lookup.  The travel-time table
    is stored in a sqlalchemy data table
    d_km=Column(Float) # Distance in kilometers
    delta=Column(Float) # Distance in degrees
    p_tt=Column(Float) # P travel-time
    s_tt=Column(Float) # S travel-time
    s_p=Column(Float) # S-P time
    """
    self.d_km=d_km
    self.delta=delta
    self.p_tt=p_tt
    self.s_tt=s_tt
    self.s_p=s_p
    
  def __repr__(self):
    return "<TTtable (%.2f %.2f %.2f %.2f %.2f)>" % (self.d_km,self.delta,self.p_tt,self.s_tt,self.s_p)


class Station1D(BaseTT1D):
  __tablename__="stations"
  id=Column(Integer,primary_key=True)
  sta=Column(String(5))
  net=Column(String(2))
  loc=Column(String(2))
  latitude=Column(Float)
  longitude=Column(Float)
  elevation=Column(Float)
  starttime=Column(DateTime)
  endtime=Column(DateTime)
  
  def __init__(self,sta,net,loc,latitude,longitude,elevation):
    self.sta=sta
    self.net=net
    self.loc=loc
    self.latitude=latitude
    self.longitude=longitude
    self.elevation=elevation
    self.starttime=None
    self.endtime=None
  
  def __repr__(self):
    return "Station <%s.%s.%s.%s %s %s %s>" % (self.sta,self.net,self.loc,self.latitude,self.longitude,self.starttime,self.endtime)    
  