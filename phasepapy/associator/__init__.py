""" phasepapy.associator

This python package contains modules for associating 
earthquake phase arrivals.

https://github.com/austinholland/PhasePApy

assoc1d
tables1d
tt_stations_1D
assoc3D
tables3D
tt_stations_3D
"""
__all__=['assoc1D','assoc3d','tables1D','tables3D','tt_stations_1d','tt_stations_3D','plot1D','plot3D']

from . import assoc1D 
#import assoc3D
#from assoc1D import *
#from assoc3D import *
#from func1D import *
#from func3D import *
from . import tables1D 
#from tables3D import *
#from search import *
from . import tt_stations_1D
from . import plot1D
#from tt_stations_3D import *