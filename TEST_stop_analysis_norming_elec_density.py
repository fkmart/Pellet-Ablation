import numpy as np 
import scipy.integrate as spint 
from gen_var import dr,pel_pot, rp, rc
import scipy.interpolate as spit
import os 
import romberg as ro 

"""This function does not normalise the density that has accummulated from all the charge
only the most recent flux of electrons to enter the region. This will need to be adapted and
modified to make life easier for myself at a later time."""
def renorm_dens(r,faux_dens,e_bins, thing3,i,r_full):
    """
    ind = np.where(thing3[:,1]< rp[i])
    x = ind[0]
    x = x[0]
    sum1 = np.sum(e_bins[x:]) # this is the fraction of particles STRIKING the pellet"""
    sum1 = 1.0 - np.sum(faux_dens) # this is the fraction of particles STRIKING the pellet
    frac_left = 1.0 - sum1

    #need to interpolate to all interior points
    rpi = rp[i]
    rci = rc[i]
    f = spit.interp1d(r, faux_dens,kind = 'cubic', fill_value = 'extrapolate')
    g = spit.PchipInterpolator(np.flip(r, axis = 0), np.flip(faux_dens, axis = 0), extrapolate = 'true')
    low_int = next(p[0] for p in enumerate(r_full) if p[1] > r[-1]) - 1
    up_int = next(p[0] for p in enumerate(r_full) if p[1] > r[0])

    low = next(p[0] for p in enumerate(r_full) if p[1] > rpi)
    up = next(p[0] for p in enumerate(r_full) if p[1] > rci)

    r1 = r_full[low_int]
    r2 = r_full[up_int]

    n_romb = len(r_full)
    """i = 0
    while (n_romb < len(r)):
        i +=1 
        n_romb = 2**i + 1"""
    ##################################################################
    int_grid = np.linspace(r1,r2, endpoint = 'true', num = n_romb)
    int_grid = np.arange(r1,r2, step = r_full[1])
    faux_dens_part = f(r_full[low_int:up_int])
    faux_dens_part = g(r_full[low_int:up_int])
    faux_dens_full = np.zeros(n_romb)
    faux_dens_full[low_int:up_int] = faux_dens_part[:]
######################################################################
    integrated = spint.romb(faux_dens_full, dx = -r[1] + r[0])
    integrated_own = ro.romberg_samp(faux_dens_full,r_full)
    real_dens = faux_dens[:]/integrated
    real_dens_full = faux_dens_full/integrated
    check = spint.romb(real_dens_full, dx = -r[1] + r[0])
    real_dens *= frac_left
    real_dens = np.append(real_dens,r)
    real_dens = np.reshape(real_dens, (int(len(real_dens)/2),2),order = 'F')
    #rdf = np.append(real_dens_full*frac_left, int_grid)
    rdf = np.append(real_dens_full*frac_left, r_full)
    rdf = np.reshape(rdf, (int(len(rdf)*0.5),2), order = 'F')
    fd_out = np.append(faux_dens_part, r_full[low_int:up_int])
    fd_out = np.reshape(fd_out, (int(0.5*len(fd_out)),2), order = 'F')
    return real_dens, integrated, rdf, fd_out