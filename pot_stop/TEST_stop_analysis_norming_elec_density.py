import numpy as np 
import scipy.integrate as spint 
from gen_var import dr,pel_pot, rp, rc
import scipy.interpolate as spit
import os 

"""This function does not normalise the density that has accummulated from all the charge
only the most recent flux of electrons to enter the region. This will need to be adapted and
modified to make life easier for myself at a later time."""
def renorm_dens(r,faux_dens,e_bins, thing3,i,r_full):
    """
    ind = np.where(thing3[:,1]< rp[i])
    x = ind[0]
    x = x[0]
    sum1 = np.sum(e_bins[x:]) # this is the fraction of particles STRIKING the pellet"""
    faux_dens = faux_dens[:-1]
    r = r[:-1]
    sum1 = 1.0 - np.sum(faux_dens) # this is the fraction of particles STRIKING the pellet
    frac_left = 1.0 - sum1

    #need to interpolate to all interior points
    rpi = rp[i]
    rci = rc[i]
    f = spit.interp1d(r[:-1], faux_dens[:-1],kind = 'cubic', fill_value = 'extrapolate') # cubic interpolation
    g = spit.PchipInterpolator(np.flip(r, axis = 0), np.flip(faux_dens, axis = 0), extrapolate = 'true') #pchip interpolation
    low = next(p[0] for p in enumerate(r_full) if p[1] > r[-1]) #index of innermost point in gen r
    up = next(p[0] for p in enumerate(r_full) if p[1] > r[0]) # index of outermost point in gen r
    r1 = r_full[low] # value of those points
    r2 = r_full[up]

    n_romb = 0
    i = 0
    while (n_romb < len(r)): # defining a romberg grid length with more points than in the stopped file
        i +=1 
        n_romb = 2**i + 1
    
    int_grid = np.linspace(r1,r2, endpoint = 'true', num = n_romb) # define the romberg grid
    faux_dens_part = f(int_grid)
    faux_dens_part = g(int_grid) # both lines interpolate using respective functions
    faux_dens_full = np.zeros(n_romb) # slightly redundant code - effectively renames  interpolated bins
    faux_dens_full[:] = faux_dens_part[:]

    integrated = spint.romb(faux_dens_full, dx = -r[1] + r[0]) # romberg integrates across the grid
    real_dens = faux_dens[:]/integrated #normalises bins to 1 by dividing through by result
    real_dens_full = faux_dens_full/integrated
    check = spint.romb(real_dens_full, dx = -r[1] + r[0]) # check 1 results to sufficient accuracy
    real_dens *= frac_left #multiply by fraction of EEDF remaining in cloud
    real_dens = np.append(real_dens,r)
    real_dens = np.reshape(real_dens, (int(len(real_dens)/2),2),order = 'F')
    rdf = np.append(real_dens_full*frac_left, int_grid)
    rdf = np.reshape(rdf, (int(len(rdf)*0.5),2), order = 'F')
    return real_dens, integrated, rdf