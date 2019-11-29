import numpy as np 
import scipy.integrate as spint 
from gen_var import dr,pel_pot, rp, rc
import scipy.interpolate as spit
import os


"""This function does not normalise the density that has accummulated from all the charge
only the most recent flux of electrons to enter the region. This will need to be adapted and
modified to make life easier for myself at a later time."""
def renorm_dens(r,faux_dens,e_bins, thing3,i,r_full):
    """ind = np.where(thing3[:,1]> rp[i])
    x = ind[-1]
    x = x[-1]
    sum1 = np.sum(faux_dens[x+1:]) # this is the fraction of particles STRIKING the pellet
    frac_left = 1.0 - sum1"""
    sum1 = 1.0 - np.sum(faux_dens) # this is the fraction of particles striking the pellet
    frac_left = 1.0 - sum1

    #need to interpolate to all interior points
    rpi = rp[i]
    rci = rc[i]
    f = spit.interp1d(r, faux_dens,kind = 'quadratic')
    
    low = next(p[0] for p in enumerate(r_full) if p[1] > r[-1])
    up = next(p[0] for p in enumerate(r_full) if p[1] > r[0])

    faux_dens_part = f(r_full[low:up])
    faux_dens_full = np.zeros(len(r_full))
    faux_dens_full[low:up] = faux_dens_part[:]

    #integrated = spint.simps(faux_dens,r,dx = dr)
    integrated = spint.simps(faux_dens_full, dx = -r[1] + r[0])
    real_dens = faux_dens[:]/integrated
    real_dens_full = faux_dens_full/integrated
    check = spint.simps(real_dens_full, dx = -r[1] + r[0])
    real_dens *= frac_left
    real_dens = np.append(real_dens,r)
    real_dens = np.reshape(real_dens, (int(len(real_dens)/2),2),order = 'F')
    return real_dens
