from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base


Base=declarative_base()

class Pick(Base):
  __tablename__="picks"
  id=Column(Integer,primary_key=True)
  sta=Column(String(5))
  chan=Column(String(3))
  net=Column(String(2))
  loc=Column(String(2))
  time=Column(DateTime)
  snr=Column(Float)
  phase=Column(String(1))
  uncert=Column(Float)
  polarity=Column(String(1))
  locate_flag=Column(Boolean)
  assoc_id=Column(Integer)
  modified_id=Column(Integer)
  t_create=Column(DateTime)
  
  def __init__(self,scnl,picktime,polarity,snr,uncert,t_create):
    self.sta=scnl.station
    self.chan=scnl.channel
    self.net=scnl.network
    self.loc=scnl.location
    self.time=picktime
    self.snr=snr
    self.uncert=uncert
    self.modified_id=None
    self.phase=None
    self.polarity=polarity
    self.locate_flag=None
    self.assoc_id=None
    self.t_create=t_create
    
  def __repr__(self):
    return "Pick <%s.%s.%s.%s %s %s %s %s>" % (self.sta,self.chan,self.net,self.loc,self.time.isoformat("T"),self.phase,self.modified_id,self.assoc_id)

class PickModified(Base):
  __tablename__="picks_modified"
  id=Column(Integer,primary_key=True)
  sta=Column(String(5))
  chan=Column(String(3))
  net=Column(String(2))
  loc=Column(String(2))
  time=Column(DateTime)
  phase=Column(String(1))
  error=Column(Float)
  locate_flag=Column(Boolean)
  assoc_id=Column(Integer)
  
  def __init__(self,sta,chan,net,loc,picktime,phase,uncert,assoc_id):
    self.sta=sta
    self.chan=chan
    self.net=net
    self.loc=loc
    self.time=picktime
    self.phase=phase
    self.uncert=uncert
    self.locate_flag=None
    self.assoc_id=assoc_id
    
  def __repr__(self):
    return "Pick <%s.%s.%s.%s %s %s %s>" % (self.sta,self.chan,self.net,self.loc,self.time.isoformat("T"),self.phase,self.assoc_id)

class Candidate(Base):
  __tablename__="candidate"
  id=Column(Integer,primary_key=True)
  ot=Column(DateTime)
  sta=Column(String(5))
  d_km=Column(Float)
  delta=Column(Float)
  weight=Column(Float)
  # P and S travel times are not completely necessary because they can be calculated but simpler to save
  tp=Column(DateTime)
  p_modified_id=Column(Integer)  # modified pick ID
  ts=Column(DateTime)
  s_modified_id=Column(Integer)  # modified pick ID
  locate_flag=Column(Boolean)
  assoc_id=Column(Integer)
  
  def __init__(self,ot,sta,d_km,delta,tp,p_modified_id,ts,s_modified_id):
    self.ot=ot
    self.sta=sta
    self.d_km=d_km
    self.delta=delta
    self.weight=None
    self.tp=tp
    self.ts=ts
    self.p_modified_id=p_modified_id
    self.s_modified_id=s_modified_id
    self.locate_flag=None
    self.assoc_id=None
  
  def __repr__(self):
    return "Candidate Event <%s %s %.2f %.2f %d %d>" % (self.ot.isoformat("T"),self.sta,self.d_km,self.delta, self.p_modified_id, self.s_modified_id)

  def __str__(self):
    return "Candidate Event <%s %s %.2f %.2f %d %d>" % (self.ot.isoformat("T"),self.sta,self.d_km,self.delta, self.p_modified_id, self.s_modified_id)

  def set_assoc_id(self,assoc_id,session,FT):
    self.assoc_id=assoc_id
    self.locate_flag=FT
    # Assign phases to modified picks
    
    # Actually only one pick_p and pick_s
    pick_p=session.query(PickModified).filter(PickModified.id==self.p_modified_id)
    for pick in pick_p:
      pick.phase='P'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT
    
    pick_s=session.query(PickModified).filter(PickModified.id==self.s_modified_id)
    for pick in pick_s:
      pick.phase='S'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT
    
    # Assign the phases to picks contribute to a modified picks
    picks_p=session.query(Pick).filter(Pick.modified_id==self.p_modified_id).all()
    for pick in picks_p:
      pick.phase='P'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT
    
    picks_s=session.query(Pick).filter(Pick.modified_id==self.s_modified_id).all()
    for pick in picks_s:
      pick.phase='S'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT  

    
class Associated(Base):
  __tablename__="associated"
  id=Column(Integer,primary_key=True)
  ot=Column(DateTime)
  ot_uncert=Column(Float)
  latitude=Column(Float)
  longitude=Column(Float)
  loc_uncert=Column(Float)
  nsta=Column(Integer)
  t_create=Column(DateTime)
  t_update=Column(DateTime)                          
  
  def __init__(self,ot,ot_uncert,latitude,longitude,loc_uncert,nsta,t_create,t_update):
    self.ot=ot
    self.ot_uncert=ot_uncert
    self.latitude=latitude
    self.longitude=longitude
    self.loc_uncert=loc_uncert
    self.nsta=nsta
    self.t_create=t_create
    self.t_update=t_update
    
  def __repr__(self):
    return "Associated Event <%s %s %.3f %.3f %d>" % (self.ot.isoformat("T"),self.ot,self.latitude,self.longitude,self.nsta)
