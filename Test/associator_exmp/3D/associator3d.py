#!/usr/bin/env python

from PhasePApy.PhasePicker import *
from obspy.core import read
import glob
import time

from sqlalchemy.orm import *
from sqlalchemy import create_engine
from PhasePApy.Associator import tables3D, assoc3D, plot3D
from datetime import datetime

# databases
db_assoc = 'sqlite:///associator3D.db'
db_tt = 'sqlite:///tt_stations_3D.db'

# connecting
engine_assoc = create_engine(db_assoc, echo = False)
tables3D.Base.metadata.create_all(engine_assoc) 
Session = sessionmaker(bind = engine_assoc)
session = Session()

# glob data
files = glob.glob('../../20130616153750/*.msd')

for file in files:
  st = read(file)
  for i in range(len(st)):
    tr = st[i]; # print tr.stats
    #=======================================================================
    # FBpicker example
    tr.detrend('linear')
    now = time.time()
    chenPicker = FBPicker(t_long = 5, freqmin = 1, mode = 'rms', t_ma = 20, nsigma = 6, \
    t_up = 0.78, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)
    scnl, picks, polarity, snr, uncert = chenPicker.picks(tr)
    print "It took %fs to pick." %(time.time()-now)
    
    # populate database with picks
    if not picks: continue
    t_create = datetime.utcnow()
    for i in range(len(picks)):
        pick = picks[i].datetime
        new_line = tables3D.Pick(scnl, pick, polarity[i], snr[i], uncert[i], t_create)
        session.add(new_line)
        session.commit()       

# associate picks with phase types
chenAssoc = assoc3D.LocalAssociator(db_assoc, db_tt, max_km = 350, aggregation = 1, aggr_norm = 'L2', assoc_ot_uncert = 3, nsta_declare = 3, nt = 31, np = 41, nr = 5)
chenAssoc.id_candidate_events()
chenAssoc.associate_candidates()

# Plot
plt = plot3D.Plot(db_assoc, db_tt)
plt.cluster_plot(assoc_ot_uncert = 3)
plt.event_plot(1)
plt.section_plot(1, files)
