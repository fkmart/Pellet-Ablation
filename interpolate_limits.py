import numpy as np
import scipy.interpolate as spint 

from gen_var import dr    

def interpol(dens_file):
    dens = dens_file[:,0] 
    space = dens_file[:,1] 

    f = spint.interp1d(space,dens, kind = 'cubic')
    low = np.amin(space)
    high = np.amax(space)
    new_space = np.arange(low, high+dr,dr)
    new_dens = f(new_space)
    return new_dens, new_space

