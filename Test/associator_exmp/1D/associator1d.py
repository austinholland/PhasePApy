#!/usr/bin/env python

from PhasePApy.PhasePicker import *
from obspy.core import read
import glob
import time

from sqlalchemy.orm import *
from sqlalchemy import create_engine
from PhasePApy.Associator import tables1D, assoc1D, plot1D
from datetime import datetime

# databases
db_assoc = 'sqlite:///associator1D.db'
db_tt = 'sqlite:///tt_stations_1D.db'

# connecting
engine_assoc = create_engine(db_assoc, echo = False)
tables1D.Base.metadata.create_all(engine_assoc) 
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
        new_line = tables1D.Pick(scnl, pick, polarity[i], snr[i], uncert[i], t_create)
        session.add(new_line)
        session.commit()       

# associate picks with phase types
chensAssoc = assoc1D.LocalAssociator(db_assoc, db_tt, max_km = 350, aggregation = 1, aggr_norm = 'L2', cutoff_outlier = 10, assoc_ot_uncert = 7, nsta_declare = 4, loc_uncert_thresh = 0.2)
candidate = chensAssoc.id_candidate_events()
chensAssoc.associate_candidates()
chensAssoc.single_phase()

# Plot
plt = plot1D.Plot(db_assoc, db_tt)
plt.cluster_plot(assoc_ot_uncert = 3)
plt.event_plot(1)
plt.section_plot(1, files)