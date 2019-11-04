import numpy as np
import scipy.integrate as spinteg
import matplotlib.pyplot as plt
import scipy.interpolate as spint
import number
import os
from electron import RME, M_fac
from MB_calc import MB
from gen_var import rc, rp, dr, e , r0
import TEST


def den_test(t):
    particle = 'electron'
    rcloud = rc[t]
    rpellet = rp[t]
    file = np.loadtxt(os.path.join('Analysed Outputs','GNU_'+particle+'_density_cloud_1keV_fs_t' +str(t) +'.txt'))
    den = file[:,1]
    check1 = np.sum(den)
    rel_tol = 0.000001
    r = np.arange(0,rcloud-rpellet,dr)
    l = len(r)
  

    """Define Maxwellian and where particles should stop, yielding a density profile that is a reflected Maxwellian"""

    den_new = spint.pchip_interpolate(np.flip(file[1:,0],axis = 0), np.flip(file[1:,1],axis = 0), r)

    """Now need number of charges that have passed through cylinder face to enter cloud"""

    no = number.number(0)

    """CHECK - to get fraction of particles stopped in the cloud from the background do the following:"""

    z = np.arange(100.0, 20000.0, 100)
    MB_dist = MB(z, 1000.0)

    sum_check = np.sum(den)
    sumthing = np.zeros(len(MB_dist))
    for i in range(0, len(MB_dist)-1):
        sumthing[i] = np.sum(MB_dist[:i])
    diff = sumthing - sum_check
    ind = np.abs(diff).argmin()

    stop_frac = spinteg.simps(MB_dist[:ind],z[:ind])
    #charge density then follows from that

    q_dens = np.copy(den_new[1:])

    """IF USING IONS THEN USE THESE FEW LINES
    depth = 5
    q_dens = np.zeros(len(den_new[1:]))
    q_dens[-depth:] = 1.0
    tot = spinteg.simps(q_dens[-depth:], r[-depth:])
    q_dens /= 1.0*tot #/(e*10**20) 
    """

    q_dens = den_new[1:].copy()
    tot = spinteg.simps(q_dens, r[1:])
    q_dens /= -1.0*tot

    """Reference values"""
    phi_0 = RME*M_fac
    epsilon_0 = 8.854*10**(-12)
    norm = e*no*(r0)**2/(epsilon_0*phi_0) # final normalisation constant

    """Now proceed with the integration"""
    phi = np.zeros(l-2)


    """Test case values to confirm validity of solver"""

    """r*= 1.0
    #q_dens = np.ones(l-1)
    q_dens /= spinteg.simps(q_dens, r[1:])
    print(spinteg.simps(q_dens, r[1:]))"""
    return q_dens, r[1:], norm

den1 , r1, norm = den_test(1)
den2, r2, norm =  den_test(2)
den3, r3, norm = den_test(3)

den_neut1 = 0.01*((1.0 + rp[1]**2)/(1.0 + (r1+rp[1])**2))
den_neut2 = 0.01*((1.0 + rp[2]**2)/(1.0 + (r2 + rp[2])**2))


def phi_calc(t, q_dens):
    rcloud = rc[t]
    rpellet = rp[t]
    r = np.arange(0,rcloud-rpellet,dr)
    phi_0 = RME*M_fac
    bc1 = 0.0
    bc2 = 0.0
    bc1 /=phi_0
    bc2 /= phi_0
    l = len(r)
    l = len(r)
    phi = np.zeros(l-2)
    A = np.zeros((l-2,l-2))
            # setting up discretisation matrix
    for i in range(0,l-2):
        A[i,i] = -2.0
        #this could be tidier but not essential to change
    for i in range(1,l-2):
        A[i,i-1] = 1.0     
    for i in range(0,l-3):
        A[i,i+1] = 1.0  
    phi = TEST.new(A, phi, q_dens, bc1, bc2)
    return phi



    
