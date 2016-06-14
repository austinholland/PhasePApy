class SCNL():
  """ This class helps us work with SCNL definitions and compare and go between different
  representations"""
  def __init__(self,input=None):
#    print input,type(input)
    if not isinstance(input, SCNL):
      self.station=None
      self.channel=None
      self.network=None
      self.location=None
    if type(input) is str:
      self.parse_scnlstr(input)
    if type(input) is list:
      if len(input)==4:
        self.station,self.channel,self.network,self.location=input
      if len(input)==3:
        self.station,self.channel,self.network=input
      if len(input)<3:
        raise SCNL_InputError("List input has %d fields minimum of 3 required" % (len(input)))
    
  def __repr__(self):
    return self.__str__()
    
  def __str__(self):
    return "%s.%s.%s.%s" % (self.station,self.channel,self.network,self.location)
  
  def to_winston(self):
    if self.location==None or self.location=='--':
      return "%s$%s$%s" % (self.station,self.channel,self.network)
    else:
      return "%s$%s$%s$%s" % (self.station,self.channel,self.network,self.location)
      
  def to_ewscnl(self):
    if self.location==None or self.location=='--':
      return "%s.%s.%s.--" % (self.station,self.channel,self.network)
    else:
      return "%s.%s.%s.%s" % (self.station,self.channel,self.network,self.location)
      
  def to_seed(self):
    if self.location==None or self.location=='--':
      return "%s.%s.%s.." % (self.station,self.channel,self.network)
    else:
      return "%s.%s.%s.%s." % (self.station,self.channel,self.network,self.location) 
           
  def parse_scnlstr(self,scnl_str):
    if re.search('\.',scnl_str):
      # Looks like an earthworm delimited scnl
      self.from_ew(scnl_str)
    if re.search('\$',scnl_str):
      # Looks like a winston scnl
      self.from_winston(scnl_str)

  def from_ew(self,scnl_str):
    scnl=scnl_str.split('.')
    self.station=scnl[0]
    self.channel=scnl[1]
    self.network=scnl[2]
    self.location=scnl[3]
  
  def from_winston(self,scnl_str):
    scnl=scnl_str.split('$')
    self.station=scnl[0]
    self.channel=scnl[1]
    self.network=scnl[2]
    if len(scnl)==4:
      self.location=scnl[3]
    else:
      self.location=None
      
    def __eq__(self,y):
      if type(y) is str:
        tmp=y
        y=SCNL()
        y.parse_scnlstr(tmp)
      if isinstance(y,SCNL):     
        if self.to_ewscnl==y.to_ewscnl:
          return True
        else:
          return False
      else: # We were given something we can't compare
        return False

    def __neq__(self,y):
      if type(y) is str:
        tmp=y
        y=SCNL()
        y.parse_scnlstr(tmp)
      if isinstance(y,SCNL):     
        if self.to_ewscnl!=y.to_ewscnl:
          return True
        else:
          return False
      else: # We were given something we can't compare
        return False
