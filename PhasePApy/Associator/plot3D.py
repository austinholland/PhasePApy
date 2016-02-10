from obspy.core.stream import read, Stream
from obspy.core.util import gps2DistAzimuth
from obspy.core import UTCDateTime
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from sqlalchemy.orm import *
from sqlalchemy import create_engine
from tables3D import *
from tt_stations_3D import *
from datetime import *

def add_subplot_axes(ax,rect,axisbg='w'):
  fig = plt.gcf()
  box = ax.get_position()
  width = box.width
  height = box.height
  inax_position  = ax.transAxes.transform(rect[0:2])
  transFigure = fig.transFigure.inverted()
  infig_position = transFigure.transform(inax_position)    
  x = infig_position[0]
  y = infig_position[1]
  width *= rect[2]
  height *= rect[2]
  subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
  x_labelsize = subax.get_xticklabels()[0].get_size()
  y_labelsize = subax.get_yticklabels()[0].get_size()
  x_labelsize *= rect[2]**0.5
  y_labelsize *= rect[3]**0.5
  subax.xaxis.set_tick_params(labelsize=x_labelsize)
  subax.yaxis.set_tick_params(labelsize=y_labelsize)
  return subax


class Plot():

  def __init__(self, db_assoc, db_tt):
  
    # Define travel time and associator database 
    engine_assoc = create_engine(db_assoc, echo = False) 
    engine_tt_stations = create_engine(db_tt, echo = False) # create a configuration file including paths
    Session1 = sessionmaker(bind = engine_assoc) # events table 
    Session2 = sessionmaker(bind = engine_tt_stations) # traveltime table
    self.assoc_db = Session1() 
    self.tt_stations_db_3D = Session2()
    
  def cluster_plot(self, assoc_ot_uncert = 3):
#        |    |                     /\
#        |    |                    /  \          /\
#        |    | /\      /\        /    \        /  \      /\
#   _____|/\__|/  \    /  \      /      \      /    \    /  \  /\____
#        |    |    \  /    \    /        \    /      \  /    \/
#        |    |     \/      \  /          \  /        \/
#        |    |              \/            \/  

    matplotlib.rcParams["axes.labelsize"]="large"
    matplotlib.rcParams["axes.linewidth"]=2.0
    matplotlib.rcParams["xtick.major.size"]=8
    matplotlib.rcParams["ytick.major.size"]=8
    matplotlib.rcParams["ytick.minor.size"]=5
    matplotlib.rcParams["xtick.labelsize"]="large"
    matplotlib.rcParams["ytick.labelsize"]="large"
    
    dt_ot=timedelta(seconds=assoc_ot_uncert)
    candidate_ots=self.assoc_db.query(Candidate).order_by(Candidate.ot).all()
    L_ots=len(candidate_ots)
    Array=[]
    for i in range(L_ots):
      cluster=self.assoc_db.query(Candidate).filter(Candidate.ot>=candidate_ots[i].ot).filter(Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      cluster_sta=self.assoc_db.query(Candidate.sta).filter(Candidate.ot>=candidate_ots[i].ot).filter(Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      l_cluster=len(set(cluster_sta))
      Array.append((i,candidate_ots[i].ot,l_cluster,len(cluster)))

    x1=np.array(Array)[:,0]
    x2=np.array(Array)[:,1]
    y1=np.array(Array)[:,2]
    y2=np.array(Array)[:,3]
    fig=plt.figure(figsize=(15,6))
    #p1=plt.plot(x2,y1,'ro-')
    p1=plt.plot(x2,y1,'o-',c='gray')
    plt.ylim(0,max(y1)+1); plt.xlim(min(x2)-timedelta(seconds=60),max(x2)+timedelta(seconds=60))
    plt.xlabel('Time (s)',fontsize=20)
    plt.ylabel('Count',fontsize=20)
    legend=plt.legend(p1,['count'])
    plt.tight_layout()
    plt.title('Unique Candidates Cluster Analysis',fontsize=20)
    fig=plt.figure(figsize=(15,6))
    #p2=plt.plot(x2,y2,'bo-')
    p2=plt.plot(x2,y2,'o-',c='gray')
    plt.ylim(0,max(y2)+1); plt.xlim(min(x2)-timedelta(seconds=60),max(x2)+timedelta(seconds=60))
    plt.xlabel('Time (s)',fontsize=20)
    plt.legend(p2,['count'])
    plt.ylabel('Count',fontsize=20)
    plt.tight_layout()
    plt.title('Total Candidates Cluster Analysis',fontsize=20)
    plt.show()    
    
  
  def event_plot(self, assoc_id, West = -104.5, East= -94, South = 33.5, North = 37.5, deltalon = 1.0, deltalat = 1.0):
    """ Plot all the circles, stations, location and residual distribution on one map by calling the event number after the event been associated.
    """
    
    fig=plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)
    #=============================================================
    ### Map bound
    # 
    #                 North
    #
    #       West                East
    #
    #                 South
    #
    #=============================================================
    ### Basemap Module
    m = Basemap(projection='merc',llcrnrlat=South,urcrnrlat=North,llcrnrlon=West,urcrnrlon=East,lat_ts=0,resolution='c')
    m.drawcountries()
    m.fillcontinents(color='white',lake_color='aqua',zorder=1)
    #=============================================================
    # draw parallels, meridians and structures.

    #m.drawparallels(np.arange(South,North,deltalat),labels=[1,0,0,0],color='gray',dashes=[1,1e-5],labelstyle='+/-',linewidth=0.1)
    #m.drawmeridians(np.arange(West,East,deltalon),labels=[0,0,0,1],color='gray',dashes=[1,1e-5],labelstyle='+/-',linewidth=0.1)
    m.drawmapboundary(fill_color='blue') 
    m.readshapefile('/Users/chenchen/Desktop/Research/gis/Geography/county','ok_counties',drawbounds=True,linewidth=1.0,color='black',zorder=2)

    m.readshapefile('/Users/chenchen/Desktop/Research/gis/Geology/OF5-95/Faults/surface','surface',linewidth=1.25,color='gray',zorder=3)
    m.readshapefile('/Users/chenchen/Desktop/Research/gis/Geology/OF5-95/Faults/surface_ot','surface_ot',linewidth=1.25,color='gray',zorder=3)
    m.readshapefile('/Users/chenchen/Desktop/Research/gis/Geology/OF5-95/Faults/subsurface_ot','subsurface_ot',linewidth=1.25,color='gray',zorder=3)
    m.readshapefile('/Users/chenchen/Desktop/Research/gis/Geology/OF5-95/Faults/subsurface','subsurface',linewidth=1.25,color='gray',zorder=3)
    
    #=============================================================
    # plot matches and mismatches circles    
    matches=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==assoc_id).filter(Candidate.locate_flag==True).all()
    mismatches=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==assoc_id).filter(Candidate.locate_flag==False).all()
    lon_eve,lat_eve=self.assoc_db.query(Associated.longitude,Associated.latitude).filter(Associated.id==assoc_id).first()
    radius_rainbow=[];s_p_rainbow=[]
    rainbow=[];sta_rainbow=[]   
    for match in matches:
      lon,lat=self.tt_stations_db_3D.query(Station3D.longitude,Station3D.latitude).filter(Station3D.sta==match.sta).first()
      s_p_rainbow.append((match.ts-match.tp).total_seconds())
      radius_rainbow.append(match.d_km)
      Color=ax._get_lines.color_cycle.next()
      rainbow.append(Color)
      sta_rainbow.append(match.sta)
      LocPair,=equi(m, lon, lat, match.d_km, lw=2., color=Color)
    legend_list={"Pair for locating":LocPair}
    
    
    for mismatch in mismatches:
      lon,lat=self.tt_stations_db_3D.query(Station3D.longitude,Station3D.latitude).filter(Station3D.sta==mismatch.sta).first()
      s_p_rainbow.append((mismatch.ts-mismatch.tp).total_seconds())
      radius_rainbow.append(mismatch.d_km)
      Color=ax._get_lines.color_cycle.next()
      rainbow.append(Color)
      sta_rainbow.append(mismatch.sta)
      UnlocPair,=equi(m, lon, lat, mismatch.d_km, ls='--', lw=1., color=Color)
    try:
      mismatch
    except NameError:
      pass
    else:
      legend_list["Pair not for locating"]=UnlocPair
      
    #=============================================================
    # plot stations
#     stations = self.assoc_db.query(Pick.sta).filter(Pick.assoc_id==assoc_id).distinct().all()  
    stations=sta_rainbow
    lon_STA=[]
    lat_STA=[]
    sta_STA=[]  
    for index,sta in enumerate(stations):  # the comma is needed
      lon,lat = self.tt_stations_db_3D.query(Station3D.longitude,Station3D.latitude).filter(Station3D.sta==sta).first()
      lon_STA.append(lon)
      lat_STA.append(lat)
      sta_STA.append(sta)
      x,y = m(lon_STA,lat_STA)
      m.scatter(x,y,marker='^',s=100,c=rainbow,zorder=3)
    for xi,yi,stai in zip(x,y,sta_STA):
      plt.text(xi,yi+1,stai,fontweight='bold',ha='center',color='k')  

    #=============================================================
    # plot event location and uncertainty
    event=self.assoc_db.query(Associated).filter(Associated.id==assoc_id).first()
    t_create=event.t_create
    t_update=event.t_update
    ot=event.ot
    lon=event.longitude
    lat=event.latitude
    location=(lon,lat)
    x,y = m(lon,lat)
    m.scatter(x,y,marker='*',s=100,c='k',zorder=3)
    #equi(m, lon, lat, uncert, color='k', lw=2.)
    #plt.title('Event %d at Origin Time: %s'%(event.id,ot), fontsize=18)
    plt.title('Event at Origin Time: %s'%(ot), fontsize=18)
    
    # add in the create and update time text
    #x_text,y_text=m(West+(East-West)*0.2,South+(North-South)*0.05)
    #plt.text(x_text,y_text,'Create Time: %s\nUpdate Time:%s'%(t_create,t_update),horizontalalignment='center',fontsize=12,fontweight='bold')

    #=============================================================
    # plot single phases circles
    single_phases=self.assoc_db.query(PickModified).filter(PickModified.assoc_id==assoc_id).filter(PickModified.locate_flag==None).all()
    for single_phase in single_phases:
      lon,lat=self.tt_stations_db_3D.query(Station3D.longitude,Station3D.latitude).filter(Station3D.sta==single_phase.sta).first()
      x,y = m(lon,lat)
      pick_time=single_phase.time
      travel_time=(pick_time-ot).total_seconds()
      if single_phase.phase=='P':
        tt,tt_uncert=singlephase(self.tt_stations_db_3D,single_phase.phase,travel_time)
        d_km=tt.d_km
        UnlocPhase,=equi(m, lon, lat, d_km, ls=':', lw=2., c='gray', zorder=4)
        m.scatter(x,y,marker='^',s=100,c='gray',zorder=5)
        plt.text(x,y+1,single_phase.sta,fontweight='bold',ha='center')  
      if single_phase.phase=='S':
        tt,tt_uncert=singlephase(self.tt_stations_db_3D,single_phase.phase,travel_time)
        d_km=tt.d_km
        UnlocPhase,=equi(m, lon, lat, d_km, ls=':', lw=2., c='gray', zorder=4)
        m.scatter(x,y,marker='^',s=100,c='gray',zorder=5)
        plt.text(x,y+1,single_phase.sta,fontweight='bold',ha='center')  
    try:
      single_phase
    except NameError:
      pass
    else:
      legend_list["Single Phase"]=UnlocPhase
    
    # add on the legend
    #legend = plt.legend(legend_list.values(),legend_list.keys(),'upper left',frameon=0) 
    
    #=============================================================
    # plot residual distribution
    subpos=[0.05,0.28,0.35,0.6]
    subax = add_subplot_axes(ax,subpos)

    for i in range(len(radius_rainbow)):
      subax.scatter(radius_rainbow[i],s_p_rainbow[i],s=50,color=rainbow[i])
    
    #print rainbow
    #subax.plot(D_S_P,S_P,'k-',linewidth=2)
    plt.xlim([0,350])
    plt.ylim([0,40])
    plt.xlabel('Distance (km)')
    plt.ylabel('S-P (s)')
    
    plt.show()  
 

 
  def section_plot(self, assoc_id, files, seconds_ahead = 5, record_length = 100, channel = 'Z'):
    
    station=self.assoc_db.query(Candidate.sta).filter(Candidate.assoc_id==assoc_id).all()
    sta_list=[]
    for sta, in station:
      sta_list.append(str(sta))
    station_single = self.assoc_db.query(Pick.sta).filter(Pick.assoc_id==assoc_id).filter(Pick.locate_flag == None).all()
    for sta, in station_single:
      sta_list.append(str(sta))
    #print sta_list
      
    eve=self.assoc_db.query(Associated).filter(Associated.id==assoc_id).first()
    # Earthquakes' epicenter
    eq_lat = eve.latitude
    eq_lon = eve.longitude
      
    # Reading the waveforms
    ST = Stream()
    for file in files:
      st = read(file)
      ST += st


    # in case of some seismometer use channel code like BH1, BH2 or BH3, resign the channel code as:
    if channel=='E' or channel=='e':
      Chan='E1'
    elif channel=='N' or channel=='n':
      Chan='N2'
    elif channel=='Z' or channel=='z':
      Chan='Z3'
    else:
      print 'Please input component E, e, N, n, Z, or z, the default is Z'
    
    # Calculating distance from headers lat/lon
    ST_new = Stream()#;print ST
    for tr in ST:
      if tr.stats.channel[2] in Chan and tr.stats.station in sta_list:
        if tr.stats.starttime.datetime < eve.ot and tr.stats.endtime.datetime > eve.ot:
          tr.trim(UTCDateTime(eve.ot-timedelta(seconds=seconds_ahead)), UTCDateTime(eve.ot+timedelta(seconds=record_length)))
          ST_new+=tr
    #print ST_new.__str__(extended=True)
 

#     remove the traces from same station, but different samples number (trace length)
#     .X31A..B Z | 2011-01-18T01:40:30.325000Z - 2011-01-18T01:41:44.975000Z | 40.0 Hz, 2987 samples
#     .WMOK..B Z | 2011-01-18T01:40:30.325000Z - 2011-01-18T01:41:44.975000Z | 40.0 Hz, 2987 samples
#     .X31A..B Z | 2011-01-18T01:40:30.325000Z - 2011-01-18T01:42:10.325000Z | 40.0 Hz, 4001 samples
#     .WMOK..B Z | 2011-01-18T01:40:30.325000Z - 2011-01-18T01:42:10.325000Z | 40.0 Hz, 4001 samples

    while True:
      ST_new_sta=[]
      for tr in ST_new:
        ST_new_sta.append(tr.stats.station)
      duplicate=list(set([tr for tr in ST_new_sta if ST_new_sta.count(tr)>1]))
      if not duplicate:
        break
      index=[i for (i,j) in enumerate(ST_new_sta) if j==duplicate[-1]]
      i=0
      while True:
        if ST_new[index[i]].stats.npts<ST_new[index[i+1]].stats.npts:
          del ST_new[index[i]]
          break
        elif ST_new[index[i]].stats.npts>=ST_new[index[i+1]].stats.npts:
          del ST_new[index[i+1]]
          break
    #print ST_new.__str__(extended=True)     


    ST_new.detrend('demean')
#     ST_new.filter('bandpass', freqmin=0.1, freqmax=100)

    factor=10
    numRows=len(ST_new)
    segs=[];ticklocs=[];sta=[];circle_x=[];circle_y=[];segs_picks=[];ticklocs_picks=[]
    for tr in ST_new:
      dmax=tr.data.max()
      dmin=tr.data.min()
      data=tr.data/(dmax-dmin)*factor
      t=np.arange(0,round(tr.stats.npts/tr.stats.sampling_rate/tr.stats.delta))*tr.stats.delta # due to the float point arithmetic issue, can not use "t=np.arange(0,tr.stats.npts/tr.stats.sampling_rate,tr.stats.delta)"
      segs.append(np.hstack((data[:,np.newaxis],t[:,np.newaxis])))
      lon,lat = self.tt_stations_db_3D.query(Station3D.longitude,Station3D.latitude).filter(Station3D.sta==tr.stats.station).first()
      distance = int(gps2DistAzimuth(lat,lon,eq_lat,eq_lon)[0]/1000.)  #gps2DistAzimuth return in meters, convert to km by /1000
#       distance=self.assoc_db.query(Candidate.d_km).filter(Candidate.assoc_id==assoc_id).filter(Candidate.sta==tr.stats.station).first()[0]#;print distance,tr.stats.station
      ticklocs.append(distance)
      sta.append(tr.stats.station)
      # DOT plot where picks are picked, notice that for vertical trace plot p is queried from Pick table, s from PickModified table
      # horizontal trace plot p and s queried from PickModified table
      if channel=='Z3':
        picks_p=self.assoc_db.query(Pick.time).filter(Pick.assoc_id==assoc_id).filter(Pick.sta==tr.stats.station).filter(Pick.chan==tr.stats.channel).filter(Pick.phase=='P').all()
        if not picks_p:
          picks_p=self.assoc_db.query(PickModified.time).filter(PickModified.assoc_id==assoc_id).filter(PickModified.sta==tr.stats.station).filter(PickModified.phase=='P').all()
        picks_s=self.assoc_db.query(PickModified.time).filter(PickModified.assoc_id==assoc_id).filter(PickModified.sta==tr.stats.station).filter(PickModified.phase=='S').all()
#         print picks_p,picks_s
      else:
        picks_p=self.assoc_db.query(PickModified.time).filter(PickModified.assoc_id==assoc_id).filter(PickModified.sta==tr.stats.station).filter(PickModified.phase=='P').all()
        picks_s=self.assoc_db.query(PickModified.time).filter(PickModified.assoc_id==assoc_id).filter(PickModified.sta==tr.stats.station).filter(PickModified.phase=='S').all()
#         print picks_p,picks_s
      picks=picks_p+picks_s
#       picks=self.assoc_db.query(PickModified.time).filter(PickModified.assoc_id==assoc_id).filter(PickModified.sta==tr.stats.station).all()
      for pick, in picks:
        index=int((pick-eve.ot+timedelta(seconds=seconds_ahead)).total_seconds()/tr.stats.delta)#;print pick,eve.ot,index,len(data)
        circle_x.append(distance+data[index])
        circle_y.append(t[index])
        # BAR plot where picks are picked  
        t_picks=np.array([t[index],t[index]])
        data_picks=np.array([data.min(),data.max()])
        segs_picks.append(np.hstack((data_picks[:,np.newaxis],t_picks[:,np.newaxis])))
        ticklocs_picks.append(distance)
    tick_max=max(ticklocs)
    tick_min=min(ticklocs)
    offsets=np.zeros((numRows,2),dtype=float)
    offsets[:,0]=ticklocs
    offsets_picks=np.zeros((len(segs_picks),2),dtype=float)
    offsets_picks[:,0]=ticklocs_picks
    
    #lines=LineCollection(segs,offsets=offsets,transOffset=None,linewidths=.25,colors=[colorConverter.to_rgba(i) for i in ('b','g','r','c','m','y','k')]) #color='gray'
    lines=LineCollection(segs,offsets=offsets,transOffset=None,linewidths=.25,color='gray')
    #lines_picks=LineCollection(segs_picks,offsets=offsets_picks,transOffset=None,linewidths=1,color='r')
    lines_picks=LineCollection(segs_picks,offsets=offsets_picks,transOffset=None,linewidths=1,color='k')
    
    #print sta,ticklocs
    fig=plt.figure(figsize=(15,8))
    ax1 = fig.add_subplot(111)
    #ax1.plot(circle_x,circle_y,'o')  # blue dots indicating where to cross the waveforms
    ax1.plot(circle_x,circle_y,'o',c='gray')
    x0 = tick_min-(tick_max-tick_min)*0.1
    x1 = tick_max+(tick_max-tick_min)*0.1
    ylim(0,record_length);xlim(0,x1)
    ax1.add_collection(lines)
    ax1.add_collection(lines_picks)
    ax1.set_xticks(ticklocs)
    ax1.set_xticklabels(sta)
    ax1.invert_yaxis()
    ax1.xaxis.tick_top()
#     ax2 = ax1.twiny()
#     ax2.xaxis.tick_bottom()   
    plt.setp(plt.xticks()[1], rotation=45)
    #xlabel('Station (km)')
    xlabel('channel: '+channel, fontsize=18)
    ylabel('Record Length (s)', fontsize=18)
#     plt.title('Section Plot of Event at %s'%(tr.stats.starttime))
#     plt.tight_layout()
    
    plt.show()


def shoot(lon, lat, azimuth, maxdist = None, Radius_earth = 6371):
  """Shooter Function
  Original javascript on http://williams.best.vwh.net/gccalc.htm
  Translated to python by Thomas Lecocq
  """
  glat1 = lat * np.pi / 180.
  glon1 = lon * np.pi / 180.
  s = maxdist / 1.852
  faz = azimuth * np.pi / 180.

  EPS= 0.00000000005
  if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
      alert("Only N-S courses are meaningful, starting at a pole!")

#   a=6378.13/1.852
  a=Radius_earth/1.852
  f=1/298.257223563
  r = 1 - f
  tu = r * np.tan(glat1)
  sf = np.sin(faz)
  cf = np.cos(faz)
  if (cf==0):
      b=0.
  else:
      b=2. * np.arctan2 (tu, cf)

  cu = 1. / np.sqrt(1 + tu * tu)
  su = tu * cu
  sa = cu * sf
  c2a = 1 - sa * sa
  x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
  x = (x - 2.) / x
  c = 1. - x
  c = (x * x / 4. + 1.) / c
  d = (0.375 * x * x - 1.) * x
  tu = s / (r * a * c)
  y = tu
  c = y + 1
  while (np.abs (y - c) > EPS):

      sy = np.sin(y)
      cy = np.cos(y)
      cz = np.cos(b + y)
      e = 2. * cz * cz - 1.
      c = y
      x = e * cy
      y = e + e - 1.
      y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
            d / 4. - cz) * sy * d + tu

  b = cu * cy * cf - su * sy
  c = r * np.sqrt(sa * sa + b * b)
  d = su * cy + cu * sy * cf
  glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
  c = cu * cy - su * sy * cf
  x = np.arctan2(sy * sf, c)
  c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
  d = ((e * cy * c + cz) * sy * c + y) * sa
  glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    

  baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

  glon2 *= 180./np.pi
  glat2 *= 180./np.pi
  baz *= 180./np.pi

  return (glon2, glat2, baz) 

def equi(m, centerlon, centerlat, radius, *args, **kwargs):
  """ http://www.geophysique.be/2011/02/20/matplotlib-basemap-tutorial-09-drawing-circles/
  """
  glon1 = centerlon
  glat1 = centerlat
  X = []
  Y = []
  for azimuth in range(0, 360):
      glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
      X.append(glon2)
      Y.append(glat2)
  X.append(X[0])
  Y.append(Y[0])

  #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
  X,Y = m(X,Y)
  plt.plot(X,Y,**kwargs)  
  return plt.plot(X,Y,**kwargs) # just return for legend    
    