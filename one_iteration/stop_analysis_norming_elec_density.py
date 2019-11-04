import numpy as np 
import scipy.integrate as spint 
from gen_var import dr, t_static,pel_pot

import os


"""This function does not normalise the density that has accummulated from all the charge
only the most recent flux of electrons to enter the region. This will need to be adapted and
modified to make life easier for myself at a later time."""
def renorm_dens(r,faux_dens,e_bins, thing3):

    ind = np.where(thing3[:,1]==0)
    x = ind[0]
    x = x[0]
    sum1 = np.sum(e_bins[x:]) # this is the fraction of particles STRIKING the pellet
    frac_left = 1.0 - sum1

    integrated = spint.simps(faux_dens,r,dx = dr)
    real_dens = faux_dens[:]/integrated
    real_dens *= frac_left
    real_dens = np.append(real_dens,r)
    real_dens = np.reshape(real_dens, (int(len(real_dens)/2),2),order = 'F')
    return real_dens 
