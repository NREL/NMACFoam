import numpy as np
from sys import argv

#geometry ========
channel_L      = 35               # length
sp_dia         = 0.25             # spacer dia
channel_W      = 0.5               # channel width
channel_span   = 5                # spanwise length
sp_gap         = 2.25        # spacing between spacers
sp_off1        = 1.75              # initial offset
nspacers       = 15        # number of spacers

sp_R = sp_dia/2.0;
sp_R_x=sp_R/np.sqrt(2.0);
sp_R_z=sp_R/np.sqrt(2.0);

#mesh ========
nr  = int(argv[1])      
nx  = int(argv[1])
nz  = int(argv[1])
ny  = 1
