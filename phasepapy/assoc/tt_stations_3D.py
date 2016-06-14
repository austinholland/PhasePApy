from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base

BaseTT3D=declarative_base()

class TTtable3D(BaseTT3D):
  __tablename__ = 'traveltimes'
  id = Column(Integer, primary_key=True)
  sta = Column(String(5))   # station
  sgid = Column(Integer)  # source grid ID
  d_km = Column(Float)      # Distance in kilometers
  delta = Column(Float)     # Distance in degrees
  p = Column(Float)         # P travel-time
  s = Column(Float)         # S travel-time
  s_p = Column(Float)       # S-P time
  pn = Column(Float)        # Pn travel-time
  sn = Column(Float)        # Sn travel-time
  sn_pn = Column(Float)     # Sn-Pn time
  
  def __init__(self,sta, sgid, d_km, delta, p_tt, s_tt, s_p, pn_tt, sn_tt, sn_pn):
    """ TTtable(d_km,delta,p_tt,s_tt,s_p)
    Create an entry in the travel-time table for later lookup.  The travel-time table
    is stored in a sqlalchemy database
    """
    self.d_km = d_km
    self.sta = sta
    self.sgid = sgid
    self.delta = delta
    self.p = p_tt
    self.s = s_tt
    self.s_p = s_p
    self.pn = pn_tt
    self.sn = sn_tt
    self.sn_pn = sn_pn
    
  def __repr__(self):
    return "<TTtable3D (d_km%.2f delta%.2f p%.2f s%.2f s-p%.2f sn%.2f pn%.2f sn-pn%.2f)>" % (self.d_km,self.delta,self.p,self.s,self.s_p,self.pn,self.sn,self.sn_pn)


class Station3D(BaseTT3D):
  __tablename__ = "stations"
  id = Column(Integer,primary_key=True)
  sta = Column(String(5))
  net = Column(String(2))
  loc = Column(String(2))
  latitude = Column(Float)
  longitude = Column(Float)
  elevation = Column(Float)
  starttime = Column(DateTime)
  endtime = Column(DateTime)
  
  def __init__(self,sta,net,loc,latitude,longitude,elevation):
    self.sta = sta
    self.net = net
    self.loc = loc
    self.latitude = latitude
    self.longitude = longitude
    self.elevation = elevation
    self.starttime = None
    self.endtime = None
  
  def __repr__(self):
    return "<%s.%s.%s.%s %s %s %s>" % (self.sta,self.net,self.loc,self.latitude,self.longitude,self.starttime,self.endtime)

class SourceGrids(BaseTT3D):
  __tablename__ = "sourcegrids"
  id = Column(Integer,primary_key=True)
  latitude = Column(Float)
  longitude = Column(Float)
  depth = Column(Float)

  def __init__(self,latitude,longitude,depth):

    self.latitude = latitude
    self.longitude = longitude
    self.depth = depth
  
  def __repr__(self):
    return "Grid <%s, %s, %s km>" % (self.latitude,self.longitude,self.depth)
    
# Base.metadata.create_all(engine)     
