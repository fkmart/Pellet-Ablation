

import numpy as np
from gen_var import dr, rc, rp, r0, e
import os
import scipy.integrate as spinteg
import scipy.interpolate as spint
import number
from MB_calc import MB
from electron import e_dist, e_bar
import matplotlib.pyplot as plt

def e_dens(t):
    particle = 'electron'
    from electron import M_fac, RME, KE_bot, KE_top
    rcloud = rc[t]
    rpellet = rp[t]
    file = np.loadtxt(os.path.join('Analysed Outputs','GNU_'+particle+'_density_cloud_1keV_fs_t' +str(t) +'.txt'))
    den = file[:,1]
    check1 = np.sum(den)
    rel_tol = 1e-6
    r = np.arange(0,rcloud-rpellet,dr)
  
    """Define Maxwellian and where particles should stop, yielding a density profile that is a reflected Maxwellian"""

    den_new = spint.pchip_interpolate(np.flip(file[1:,0],axis = 0), np.flip(file[1:,1],axis = 0), r)
    ch = spinteg.simps(den_new, r)

    """Now need number of charges that have passed through cylinder face to enter cloud"""

    no = number.number(t)

    """CHECK - to get fraction of particles stopped in the cloud from the background do the following:"""

    #z = np.arange(100.0, 20000.0, 100)
    #MB_dist = MB(z, 1000.0)
    MB_dist = MB(e_dist, e_bar)
    """If following line returns a value of 1 then MB distribution for background is normalised"""
    checkag = spinteg.simps(MB_dist, e_dist) #""" - this is the one that is normalised, the itnegrated value"""

    stop_frac = np.sum(den)
#Stuff between dashed lines needs amended
#-------------------------------------------------------------
    
    """
    sum_check = np.sum(den)
    sumthing = np.zeros(len(MB_dist))
    for i in range(0, len(MB_dist)):
        sumthing[i] = np.sum(MB_dist[:i])
    diff = sumthing - sum_check  #- this results will ALWAYS BE WRONG MUST BE FIXED ASAP
    ind = np.abs(diff).argmin()

    stop_frac = spinteg.simps(MB_dist[:ind],e_dist[:ind])
    #charge density then follows from that"""
#-----------------------------------------------------------
    q_dens = np.copy(den_new[1:])

    q_dens = den_new[1:].copy()
    tot = spinteg.simps(q_dens, r[1:])
    tot2 = np.sum(q_dens)
    q_dens /= -1.0*tot/stop_frac
    check_tot = spinteg.simps(q_dens,r[1:])

    """Reference values"""
    phi_0 = RME*M_fac
    epsilon_0 = 8.854*10**(-12)
    norm = e*no*(r0)**2/(epsilon_0*phi_0) # final normalisation constant

    return q_dens, r[1:], norm, tot 
