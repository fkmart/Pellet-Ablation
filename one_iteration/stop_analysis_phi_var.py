import numpy as np 

from gen_var import * 

particle = 'electron'

if (particle == 'electron'):
    from electron import *
elif (particle == 'proton'):
    from proton import *
else:
    print('Particle is not known. Please create data file for particle and re-run') 

rp = stop_calc_rp_rc.calc_rp(t)
rc = stop_calc_rp_rc.calc_rc(t)

r = np.arange(rp,rc,dr)
