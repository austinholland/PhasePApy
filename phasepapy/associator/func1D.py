from .tt_stations_1D import *

def tt_km(session,d_km):
  """ tt_km(session,d_km)
  Look up travel time from tt_table for the closest distance to entered value stored in the 
  tt_table database. km_difference is the distance difference between the stored value
  returned in the lookup and the distance requested.
  
  session is the database connection
  
  Returns the TTtable_object,km_difference 
  """

  min=session.query(TTtable1D).filter(TTtable1D.d_km<=d_km).order_by(TTtable1D.d_km.desc()).first()
#   print min
  max=session.query(TTtable1D).filter(TTtable1D.d_km>=d_km).order_by(TTtable1D.d_km).first()
#   print max.d_km
  if abs(min.d_km-d_km) <= abs(max.d_km-d_km):
    return min,abs(min.d_km-d_km)
  else:
    return max,abs(max.d_km-d_km)

   
def tt_s_p(session,s_p):
  """ tt_s_p(session,s_p)
  Look up the distance for an S-P travel time observation and return the closest value
  stored in the travel-time table
  
  Returns the closest tt_object,s_p_difference
  """
#   print 's_p:',s_p
  min=session.query(TTtable1D).filter(TTtable1D.s_p<=s_p).order_by(TTtable1D.s_p.desc()).first()
#   print 'min.s_p:',min.s_p
  max=session.query(TTtable1D).filter(TTtable1D.s_p>=s_p).order_by(TTtable1D.s_p).first()
#   print 'max.s_p:',max.s_p
  if abs(min.s_p-s_p) <= abs(max.s_p-s_p):
    return min,abs(min.s_p-s_p)
  else:
    return max,abs(max.s_p-s_p)    
    
