from .tt_stations_3D import *

def tt_km(session,d_km):
  """ tt_km(session,d_km)
  Look up travel time from tt_table for the closest distance to entered value stored in the 
  tt_table database. km_difference is the distance difference between the stored value
  returned in the lookup and the distance requested.
  
  session is the database connection
  
  Returns the TTtable_object,km_difference 
  """

  min=session.query(TTtable3D).filter(TTtable3D.d_km<=d_km).order_by(TTtable3D.d_km.desc()).first()
#   print min
  max=session.query(TTtable3D).filter(TTtable3D.d_km>=d_km).order_by(TTtable3D.d_km).first()
#   print max.d_km
  if abs(min.d_km-d_km) <= abs(max.d_km-d_km):
    return min,abs(min.d_km-d_km)
  else:
    return max,abs(max.d_km-d_km)
    
    
def tt_s_p(session,s_p,uncert):
  """ tt_s_p(session,s_p,uncert)
  Look up the distance for an S-P travel time observation and return all the nodes within
  the s_p +/- uncert
  
  Returns the TTtable3D_object
  """
  
  nodes = session.query(TTtable3D).filter(TTtable3D.s_p<=s_p+uncert).filter(TTtable3D.s_p>=s_p-uncert).all()
  return nodes