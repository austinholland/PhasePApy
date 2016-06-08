from tt_stations_3D import *

def PyramidSearching(session_tt_3D, nt, np, nr, tt):
  
  session = session_tt_3D
   
  nt_seg = nt - 1
  np_seg = np - 1
  nr_seg = nr - 1
  
  # ======================================================================================
  # course search
  
  if np_seg > 3:
    np_1, np_2, np_1_cpt, np_2_cpt = Segments(np_seg) 
    np_cpt_set = [(np_1, np_1_cpt), (np_2, np_2_cpt)]
  else:
    np_cpt_set = [(np_seg, np_seg//2 + 1)]
    
  if nt_seg > 3:
    nt_1, nt_2, nt_1_cpt, nt_2_cpt = Segments(nt_seg) 
    nt_cpt_set = [(nt_1, nt_1_cpt), (nt_2, nt_2_cpt)]
  else:
    nt_cpt_set = [(nt_seg, nt_seg//2 + 1)]
    
  if nr_seg > 3:
    nr_1, nr_2, nr_1_cpt, nr_2_cpt = Segments(nr_seg) 
    nr_cpt_set = [(nr_1, nr_1_cpt), (nr_2, nr_2_cpt)]
  else:
    nr_cpt_set = [(nr_seg, nr_seg//2 + 1)]
  
    
  while nt_seg > 3 or np_seg > 3 or nr_seg > 3:
    
    #print np_cpt_set
    #print nt_cpt_set
    #print nr_cpt_set 
    
    # checking points
    #now = time.time()
    
    cpt_set = []
    
    for i in range(len(nr_cpt_set)):
      for j in range(len(nt_cpt_set)):
        for k in range(len(np_cpt_set)):
          vars()['cpt_%d'%(i * len(nt_cpt_set) * len(np_cpt_set) + j * len(np_cpt_set) + k)] \
          = (nr_cpt_set[i][1] - 1) * (np - 1) * (nt - 1) + (nt_cpt_set[j][1] - 1) * (np - 1) + (np_cpt_set[k][1] - 1)
          
          cpt_set.append((i,j,k,vars()['cpt_%d'%(i * len(nt_cpt_set) * len(np_cpt_set) + j * len(nt_cpt_set) + k)]))    
    #print 'cp_set', cpt_set
    #print 'checking points: ', time.time()-now
    
     
    # calculate rms for each checking point
    RMS = []; tt_new = []
    for i in range(len(cpt_set)):
      rms = 0
      for j in range(len(tt)):
        sta = tt[j][0]
        ttp = tt[j][1]
        tts = tt[j][2]
        ot = tt[j][3]
        # check if sta in TTtable3D
        if session.query(Station3D).filter(Station3D.sta==sta).all():
          sta_id, = session.query(Station3D.id).filter(Station3D.sta == sta).first()
          #print 'sta_id: ', sta_id
          #now = time.time() 
          id = (sta_id - 1) * (nt - 1) * (np - 1) * (nr - 1) + cpt_set[i][3] + 1
          ttp3D,tts3D = session.query(TTtable3D.p,TTtable3D.s).filter(TTtable3D.id==id).first()
          #print ttp3D,tts3D
          #print 'querying time: ', time.time()-now
          #now = time.time()
          rms += (ttp - ttp3D)**2
          rms += (tts - tts3D)**2
          if (sta,ttp,tts,ot) not in tt_new:
            tt_new.append((sta,ttp,tts,ot))
          #print 'rms calculation time: ', time.time()-now
          
      # rms of p and s together
      rms = (rms / len(tt)) ** 1/2
    
      # might worth calculating rms for each phase
      #rms_p = (rms_p / len(tt)) ** 1/2
      #rms_s = (rms_s / len(tt)) ** 1/2
    
      RMS.append(rms)
    
    #print RMS
    min_index = RMS.index(min(RMS))

    #print min_index
    
    # determine the index of the checking point set: i is indexing nr, j is indexing nr, k is indexing np
    x = cpt_set[min_index][2]
    y = cpt_set[min_index][1]
    z = cpt_set[min_index][0]  
    
    np_seg = np_cpt_set[x][0]
    np_cpt = np_cpt_set[x][1]
    if np_seg > 3:
      np_1, np_2, np_1_cpt, np_2_cpt = Segments(np_seg) 
      np_1_cpt = np_cpt - (np_1 - np_1_cpt); np_2_cpt += np_cpt - np_1
      np_cpt_set = [(np_1, np_1_cpt), (np_2, np_2_cpt)]
    else:
      np_cpt_set = [(np_seg, np_cpt)]
      X = np_cpt
      
    nt_seg = nt_cpt_set[y][0]
    nt_cpt = nt_cpt_set[y][1]
    if nt_seg > 3:
      nt_1, nt_2, nt_1_cpt, nt_2_cpt = Segments(nt_seg) 
      nt_1_cpt = nt_cpt - (nt_1 - nt_1_cpt); nt_2_cpt += nt_cpt - nt_1
      nt_cpt_set = [(nt_1, nt_1_cpt), (nt_2, nt_2_cpt)]
    else:
      nt_cpt_set = [(nt_seg, nt_cpt)]
      Y = nt_cpt
  
    nr_seg = nr_cpt_set[z][0]
    nr_cpt = nr_cpt_set[z][1]
    if nr_seg > 3:
      nr_1, nr_2, nr_1_cpt, nr_2_cpt = Segments(nr_seg) 
      nr_1_cpt = nr_cpt - (nr_1 - nr_1_cpt); nr_2_cpt += nr_cpt - nr_1
      nr_cpt_set = [(nr_1, nr_1_cpt), (nr_2, nr_2_cpt)]
    else:
      nr_cpt_set = [(nr_seg, nr_cpt)]  
      Z = nr_cpt
  
  # final grid index, not id. id = index + 1
  grid_index = cpt_set[min_index][3]
  
  #print 'X, Y, Z:', X, Y, Z
  #print 'course searching RMS: ', RMS[min_index]
  
  sourcegrid = grid_index + 1
  #print 'sourcegrid: ', sourcegrid
  
  location = session.query(SourceGrids).filter(SourceGrids.id == sourcegrid).first()
  #print 'location: ', location
  
#   if session.query(Station3D).filter(Station3D.sta==sta).all():  
#     ttp3D = session.query(TTtable3D.p).filter(TTtable3D.sta==sta).all()[grid_index][0]
#     tts3D = session.query(TTtable3D.s).filter(TTtable3D.sta==sta).all()[grid_index][0]
#     print 'ttp3D:', ttp3D, 'tts3D', tts3D
  
  # ======================================================================================
  # fine search
  
  while True: 
    # finer checking points  
    
    Y_fine = []
    if nt - 1 == 1:
      Y_fine = [Y]
    elif nt - 1 > 1:
      if Y == 1: 
        Y_fine.append(Y); Y_fine.append(Y + 1)
      elif Y == nt - 1:
        Y_fine.append(Y - 1); Y_fine.append(Y)
      else:
        Y_fine.append(Y - 1); Y_fine.append(Y); Y_fine.append(Y + 1)
    #print Y_fine
    
    X_fine = []
    if np - 1 == 1:
      X_fine = [X]
    elif np - 1 > 1:
      if X == 1: 
        X_fine.append(X); X_fine.append(X + 1)
      elif X == np - 1:
        X_fine.append(X - 1); X_fine.append(X)
      else:
        X_fine.append(X - 1); X_fine.append(X); X_fine.append(X + 1)
    #print X_fine
    
    Z_fine = []
    if nr - 1 == 1:
      Z_fine = [Z]
    elif nr - 1 > 1:
      if Z == 1: 
        Z_fine.append(Z); Z_fine.append(Z + 1)
      elif Z == nr - 1:
        Z_fine.append(Z - 1); Z_fine.append(Z)
      else:
        Z_fine.append(Z - 1); Z_fine.append(Z); Z_fine.append(Z + 1)
    #print Z_fine
    
    cpt_set_f = []
    for i in range(len(Z_fine)):
        for j in range(len(Y_fine)):
          for k in range(len(X_fine)):
            vars()['cpt_%d'%(i * len(Y_fine) * len(X_fine) + j * len(Y_fine) + k)] \
            = (Z_fine[i] - 1) * (np - 1) * (nt - 1) + (Y_fine[j] - 1) * (np - 1) + (X_fine[k] - 1)
          
            cpt_set_f.append((i,j,k,vars()['cpt_%d'%(i * len(Y_fine) * len(X_fine) + j * len(Y_fine) + k)]))    
    #print cpt_set_f
  
    # calculate rms for each checking point
    RMS_f = []
    for i in range(len(cpt_set_f)):
      rms = 0
      for j in range(len(tt)):
        sta = tt[j][0]
        ttp = tt[j][1]
        tts = tt[j][2]
        # check if sta in TTtable3D
        if session.query(Station3D).filter(Station3D.sta==sta).all():
          #now = time.time()  
          sta_id, = session.query(Station3D.id).filter(Station3D.sta == sta).first()
          id = (sta_id - 1) * (nt - 1) * (np - 1) * (nr - 1) + cpt_set_f[i][3] + 1;
          ttp3D,tts3D = session.query(TTtable3D.p,TTtable3D.s).filter(TTtable3D.id==id).first()
          #print 'querying time: ', time.time()-now
          
          #now = time.time()
          rms += (ttp - ttp3D)**2
          rms += (tts - tts3D)**2
          #print 'rms calculation time: ', time.time()-now
  
      # rms of p and s together
      rms = (rms / len(tt)) ** 1/2
  
      # might worth calculating rms for each phase
      #rms_p = (rms_p / len(tt)) ** 1/2
      #rms_s = (rms_s / len(tt)) ** 1/2
  
      RMS_f.append(rms)
  
    #print RMS_f
    min_index_f= RMS_f.index(min(RMS_f))

    #print min_index_f
    #print 'fine searching RMS: ', RMS_f[min_index_f]
    
    # determine the index of the checking point set: i is indexing nr, j is indexing nr, k is indexing np
    x_index = cpt_set_f[min_index_f][2]
    x = X_fine[x_index]
    y_index = cpt_set_f[min_index_f][1]
    y = Y_fine[y_index]
    z_index = cpt_set_f[min_index_f][0]
    z = Z_fine[z_index]
      
    #print 'x, y, z:', x, y, z
    
    if [x, y, z] == [X, Y, Z]:
      break
    else:
      [X, Y, Z] = [x, y, z]
  
  
  grid_index = cpt_set_f[min_index_f][3]
  sourcegrid = grid_index + 1
  #print 'sourcegrid: ', sourcegrid
  location = session.query(SourceGrids).filter(SourceGrids.id == sourcegrid).first()
  #print 'location: ', location
  
  return tt_new, sourcegrid, RMS_f[min_index_f]


# ========================================================================================
def Segments(n):
  """
  n has to be greater than 3
  """
  # divide Segment into 2 sub Segments
  if n % 2 == 0:
    n_1 = n / 2
    n_2 = n / 2
  else:
    n_1 = (n + 1) / 2
    n_2 = n - n_1
  
  # determine the checking point index of the two sub Segments
  if n_1 % 2 == 0:
    n_1_checking_point = n_1 / 2
  else:
    n_1_checking_point = (n_1 + 1) / 2
  
  if n_2 % 2 == 0:
    n_2_checking_point = n_2 / 2 + n_1
  else:
    n_2_checking_point = (n_2 + 1) / 2 + n_1
      
  return n_1, n_2, n_1_checking_point, n_2_checking_point
  




def SphericalSearching(nt, np, nr, tt):
  """
  Replace coarse searching with spherical searching, then do fine searching.
  """
  pass

