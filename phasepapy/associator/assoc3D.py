import numpy as np
from sqlalchemy.orm import *  # session
from sqlalchemy import create_engine

from tables3D import *
from tt_stations_3D import *
from func3D import *
from search import *
from datetime import *
from operator import itemgetter
from itertools import combinations
#import time       "from datetime import *" will import time, name space will be overwritten

class LocalAssociator():
  """
  The 3D Associator associate picks with travel time curve of 3D velocity.
  """
  def __init__(self, db_assoc, db_tt, max_km = 350, aggregation = 1, aggr_norm = 'L2', assoc_ot_uncert = 3, nsta_declare = 3, nt = 31, np = 41, nr = 5):
    """
    Parameters:
    db_assoc: associator database
    db_tt: travel time table database
    max_km: maximum distance of S-P interval in distance
    aggregation: the coefficient multiplied to minimum travel time
    aggr_norm: L2: median; L1: mean
    assoc_ot_uncert: origin time uncertainty window
    nsta_declare: minimum station number to declare a earthquake
    nt, np, nr: node geometry
    """
    
    engine_associator = create_engine(db_assoc, echo=False)
    engine_tt_stations_3D = create_engine(db_tt, echo=False)
    Base.metadata.create_all(engine_associator)               # Base is from the imported tables3D
    BaseTT3D.metadata.create_all(engine_tt_stations_3D)       # BaseTT3D is from the imported tables3D
    Session1=sessionmaker(bind = engine_associator) # events table
    Session2=sessionmaker(bind = engine_tt_stations_3D)   
    self.assoc_db                 =         Session1() 
    self.tt_stations_db_3D        =         Session2()                             
    
    self.max_km                   =         max_km
    tmp, d_diff                   =         tt_km(self.tt_stations_db_3D,self.max_km) # From max distance set our maximum travel_time
    self.max_tt                   =         tmp.s
    self.max_s_p                  =         tmp.s_p
    self.min_s_p                  =         self.tt_stations_db_3D.query(TTtable3D.s_p).order_by(TTtable3D.s_p).first()[0]
    self.aggregation              =         aggregation
    self.aggr_window              =         self.aggregation * self.min_s_p
    self.aggr_norm                =         aggr_norm         # L1 takes the mean; L2 takes the median
    self.assoc_ot_uncert          =         assoc_ot_uncert   # Uncertainty of origin times of candidates
    self.nsta_declare             =         nsta_declare      # number observation to declare an evnet
    
    self.nt                       =         nt
    self.np                       =         np
    self.nr                       =         nr
    
  def id_candidate_events(self):
    """ Create a set of possible candidate events from our picks table.
    Where session is the connection to the sqlalchemy database.
    This method simply takes all picks with time differences less than our maximum S-P
    times for each station and generates a list of candidate events.
    """
    #now1 = time.time()
    #############
    # Get all stations with unnassoiated picks
    stations=self.assoc_db.query(Pick.sta).filter(Pick.assoc_id==None).distinct().all()

    for sta, in stations:  # the comma is needed
      picks=self.assoc_db.query(Pick).filter(Pick.sta==sta).filter(Pick.assoc_id==None).order_by(Pick.time).all()
      # Condense picktimes that are within our pick uncertainty value picktimes are python datetime objects
      if stations.index((sta,))==0: #stupid tuple
        counter0=0
        picktimes_new,counter=pick_cluster(self.assoc_db,picks,self.aggr_window,self.aggr_norm,counter0)
      else:
        picktimes_new,counter=pick_cluster(self.assoc_db,picks,self.aggr_window,self.aggr_norm,counter)
      picks_modified=self.assoc_db.query(PickModified).filter(PickModified.sta==sta).filter(PickModified.assoc_id==None).order_by(PickModified.time).all()
      
      # Generate all possible candidate events
      for i in range(0, len(picks_modified) - 1):
        for j in range(i + 1,len(picks_modified)):
          s_p = (picks_modified[j].time - picks_modified[i].time).total_seconds()#; print s_p
          if s_p <= self.max_s_p and s_p >= self.min_s_p:
            nodes = tt_s_p(self.tt_stations_db_3D, s_p, .1)
            ots = []; d_kms = []; deltas = []
            for node in nodes:
              ots.append(picks_modified[i].time-timedelta(seconds=node.p))
              d_kms.append(node.d_km)
              deltas.append(node.delta)
            ot,ot_uncert = datetime_statistics(ots,norm='L2') # L1: mean, L2: median
            d_km = np.median(d_kms)
            delta = np.median(deltas)
            #print 'ot:',ot, d_km, delta
            #print 'length of nodes:',len(nodes)
    
            #ot=picks_modified[i].time-timedelta(seconds=tt.p_tt)
            new_candidate=Candidate(ot,sta,d_km,delta,picks_modified[i].time,picks_modified[i].id,picks_modified[j].time,picks_modified[j].id)
            self.assoc_db.add(new_candidate)
            self.assoc_db.commit()
    
    #print 'id_candidate time in seconds: ',time.time()-now1
      
  def associate_candidates(self):
    """ Associate all possible candidate events by comparing the projected origin-times.  At
    this point we are not dealing with the condition that more picks and candidate events 
    could be arriving while we do our event associations.
    """
    #now2 = time.time()
    
    dt_ot=timedelta(seconds=self.assoc_ot_uncert)
    
    # Query all candidate ots
    candidate_ots=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).order_by(Candidate.ot).all()
    L_ots=len(candidate_ots) #; print L_ots
    Array=[]
    for i in range(L_ots):
      cluster=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).filter(Candidate.ot>=candidate_ots[i].ot).filter(Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      cluster_sta=self.assoc_db.query(Candidate.sta).filter(Candidate.assoc_id==None).filter(Candidate.ot>=candidate_ots[i].ot).filter(Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      l_cluster=len(set(cluster_sta))
      Array.append((i,l_cluster,len(cluster)))
    #print Array
    Array.sort(key=itemgetter(1), reverse=True)  #sort Array by l_cluster, notice Array has been changed
    #print Array
    
    #print 'candidates_ots:', time.time()-now2
          
    for i in range(len(Array)):
      index=Array[i][0]
      if Array[i][1]>=self.nsta_declare:
        matches=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).filter(Candidate.ot>=candidate_ots[index].ot).filter(Candidate.ot<(candidate_ots[index].ot+dt_ot)).order_by(Candidate.ot).all() 
        
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # remove the candidates with the modified picks has been associated
        picks_associated_id=list(set(self.assoc_db.query(PickModified.id).filter(PickModified.assoc_id!=None).all()))
        index_matches=[]
        for id, in picks_associated_id:
          for i,match in enumerate(matches):
            if match.p_modified_id==id or match.s_modified_id==id:
              index_matches.append(i)        
        # delete from the end
        if index_matches:
          for j in sorted(set(index_matches),reverse=True):
            del matches[j]
        # remove the candidates with the modified picks has been associated
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        # 3D Associator
        #now = time.time()
        tt = []
        for match in matches:
          #print 'sta:',match.sta
          match_p = match.tp
          match_s = match.ts
          match_ot = match.ot
          match_ttp = (match_p - match_ot).total_seconds()
          match_tts = (match_s - match_ot).total_seconds()
          tt.append([match.sta, match_ttp, match_tts, match.ot, match])
        #print matches
        #print tt
        
        cb = self.comb(tt)
        #print 'cb',cb
        
        rms_sort = []
        for i in range(len(cb)):
          tt_cb = cb[i]
          if len(tt_cb) >= self.nsta_declare: # self.nsta_declare has to be greater than or equal to 3
            tt_new, sourcegrid, rms = PyramidSearching(self.tt_stations_db_3D, self.nt,self.np,self.nr,tt_cb)
            rms_sort.append((tt_new, sourcegrid, rms, i))
        
        # It is possible to have empty rms_sort if all the tt_cb have length less than 3
        if rms_sort:
          rms_sort.sort(key = itemgetter(2))
          tt_new, sourcegrid, rms, index = rms_sort[0]
          lat, lon, dep = self.tt_stations_db_3D.query(SourceGrids.latitude, SourceGrids.longitude, SourceGrids.depth).filter(SourceGrids.id == sourcegrid).first()
          
          nsta = len(tt_new)
          
          all_ots = []
          for i in range(nsta):  
            all_ots.append(tt_new[i][3])
          origintime, ot_unc = datetime_statistics(all_ots)
          # in 3D Associator, use rms of picks instead of loc_uncertainty
          t_create = datetime.utcnow()
          t_update = datetime.utcnow()
          new_event=Associated(origintime,round(ot_unc,3),lat,lon,dep,round(rms,3),nsta,t_create,t_update)    
          self.assoc_db.add(new_event)
          self.assoc_db.flush()
          self.assoc_db.refresh(new_event)
          self.assoc_db.commit()
          event_id=new_event.id
          print 'event_id:',event_id 
          print 'ot:', origintime, 'ot_uncert:', round(ot_unc,3), 'loc:', lat,lon, 'rms:', round(rms,3), 'nsta:', nsta
           
          for tt_tuple in cb[index]:
            match = tt_tuple[4]
            match.set_assoc_id(event_id,self.assoc_db,True)
          self.assoc_db.commit()
        #print 'time: ', time.time() - now  
        # 3D Associator
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

      else:
        break
  
  
  # create the combinations from different stations
  def comb(self,tt):
    L = len(set([item[0] for item in tt]))   # length of the set(sta)
    cb = list(combinations((tt),L))          # combinations of the array, some of them are repeated such as (sta1, sta1, sta2,...)
    
    # remove those combinations of repeated station
    index = []
    for i in range(len(cb)):
      temp = []
      for j in range(L):
        temp.append(cb[i][j][0])
      l = len(set(temp))
      if l < L:
        index.append(i) 
    index.reverse()
    for i in index:
      del cb[i]
 
    # only return combinations of different stations
    return cb   




def datetime_statistics(dt_list,norm='L2'):
  """ mean,std=datetime_statistics(datetime_list)
  Calculate the mean and standard deviations in seconds of a list of datetime values
  """
  offsets=[]
  for dt in dt_list:
    offsets.append((dt-dt_list[0]).total_seconds())
  if norm=='L1':
    mean_offsets=np.mean(offsets)
  elif norm=='L2':
    mean_offsets=np.median(offsets)
  std_offsets=np.std(offsets)
  return dt_list[0]+timedelta(seconds=mean_offsets),std_offsets




  
def pick_cluster(session,picks,pickwindow,aggr_norm,counter):
  """ cleaning up very closed picks on different channels of same station
  """
  #                     |    |                     /\
  #                     |    |                    /  \          /\
  #                     |    | /\      /\        /    \        /  \      /\
  #        _____________|/\__|/  \    /  \      /      \      /    \    /  \  /\_________
  #                     |    |    \  /    \    /        \    /      \  /    \/
  #                     |    |     \/      \  /          \  /        \/
  #                     |    |              \/            \/

  # pickwindow:          ----                                      better to set pickwindow==t_up, t_up is to clean closed picks
  # STA1 E   -----------|----|--------------------|--------------
  # STA1 N   ------------|-------------------------|-------------
  # STA1 Z   -------------|-------------------------|------------
  # stack    -----------|||--|--------------------|||------------  
  # cluster STA1 --------|---|---------------------|-------------  chen highly recommend to use norm=='L2' to lower the effect of outlier, L2 takes median
  # ARGUE: whether only take the median or mean of the picks from different stations? won't count the followings after first one
  # 
  
  picks_new=[]
  # only one pick in picks
  if len(picks)==1:
    cluster=[];cluster.append(picks[0]);cluster_time=[];cluster_time.append(picks[0].time)
    picks[0].modified_id=1+counter # assign modified id to picks
    counter+=1
    pickave,pickstd=datetime_statistics(cluster_time,aggr_norm)
    # append the row to the picks_new, not only the pick time
    picks_new.append(picks[0])
    pick_modified=PickModified(picks[0].sta,picks[0].chan,picks[0].net,picks[0].loc,picks[0].time,picks[0].phase,round(pickstd,3),picks[0].assoc_id)
    session.add(pick_modified)
    session.commit()
    
  # more than one pick in picks
  else:
    j=0
    counter=1+counter
    while True:
      i=j
      cluster=[];cluster.append(picks[i]);cluster_time=[];cluster_time.append(picks[i].time);channel=[];channel.append(picks[i].chan)
      picks[i].modified_id=counter
      while True:
        # cluster picks of different channels; notice that for the closed picks on the same station, those picks behind the first pick could be separated lonely or separated cluster
        if picks[i+1].chan not in channel and (picks[i+1].time-picks[i].time).total_seconds()<pickwindow:
          cluster.append(picks[i+1])
          cluster_time.append(picks[i+1].time)
          channel.append(picks[i+1].chan)
          picks[i+1].modified_id=counter # assign modified id to picks     
          i=i+1
          # make sure do not go over the range limit because j=i+1 below, jump out inner while loop
          if i==len(picks)-1:
            break
        # elif is dealing with the exactly same picks, probably from processing same stream twice
        elif picks[i+1].sta==picks[i].sta and picks[i+1].chan==picks[i].chan and picks[i+1].time==picks[i].time: # and picks[i+1].snr==picks[i].snr and picks[i+1].phase==picks[i].phase and picks[i+1].uncert==picks[i].uncert:
          cluster.append(picks[i+1])
          cluster_time.append(picks[i+1].time)
          channel.append(picks[i+1].chan)
          picks[i+1].modified_id=counter # assign modified id to picks     
          i=i+1
          # make sure do not go over the range limit because j=i+1 below, jump out inner while loop
          if i==len(picks)-1:
            break
        else:
          break
      pickave,pickstd=datetime_statistics(cluster_time,aggr_norm)
      
      # append whole rows to the picks_new, not only the pick time
      for pick in cluster:
        if (pick.time-pickave).total_seconds()>=0:
          break
      picks_new.append(pick)
      pick_modified=PickModified(pick.sta,pick.chan,pick.net,pick.loc,pick.time,pick.phase,round(pickstd,3),pick.assoc_id)
      session.add(pick_modified)
      session.commit()
      # next cluster
      j=i+1
      counter=counter+1
      
      # jump outer while loop and compare last two picks. For the situation that last one is ungrouped, use if statement to add in picks_new
      if j>=len(picks)-1:
        if (picks[-1].time-picks[-2].time).total_seconds()>pickwindow:
          picks_new.append(picks[-1])
          picks[-1].modified_id=counter # assign modified id to picks
          pick_modified=PickModified(picks[-1].sta,picks[-1].chan,picks[-1].net,picks[-1].loc,picks[-1].time,picks[-1].phase,round(pickstd,3),picks[-1].assoc_id)
          session.add(pick_modified)
          session.commit()
        else:
          if picks[-1] in cluster:
            counter-=1
          else:
            picks[-1].modified_id=counter
            pick_modified=PickModified(picks[-1].sta,picks[-1].chan,picks[-1].net,picks[-1].loc,picks[-1].time,picks[-1].phase,round(pickstd,3),picks[-1].assoc_id)
            session.add(pick_modified)
            session.commit()
        break     
  
  return picks_new, counter
  