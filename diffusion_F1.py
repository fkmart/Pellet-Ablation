# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 08:30:28 2020

@author: Kyle
"""

import numpy as np
import scipy.signal as scisig
import os
import romberg as ro


def diffusion(n,r,D,il,ir,dt):
    l = len(n)
    J = np.zeros(l)
    dx = r[1] - r[0]
    n[il] = 0.0
    i_max = n.argmax()
    dn = np.zeros(l)
    dd = np.zeros(l)
    for i in range(i_max, 0,-1):
        dndr_l = (n[i-1] - n[i])/dx 
        dndr_r = (n[i] - n[i+1])/dx
        J[i] = -D[i]*(dndr_l + dndr_r)
    for i in range(i_max, l-1):
        dndr_l = (n[i-1] - n[i]) / dx
        dndr_r = (n[i] - n[i+1])/dx
        J[i] = -D[i]*(dndr_l + dndr_r)
    for i in range(0, l-1): 
        dd[i] = (D[i] - D[i+1])/dx
    dn[:] = dt*J[:]
    n_out = n + dn 
    n_out[ir:] = 0.0
    n_out[:il] = 0.0
    ion_flux = np.sum(J[:il])
    return ion_flux, n_out

def diff_flux(n,r,D,il,ir,dt):
    l = len(n)
    dx = r[1] - r[0]
    dn = np.zeros(l)
    d2n = np.zeros(l)
    dd = np .zeros(l)
    term1 = np.zeros(l)
    term2 = np.zeros(l)
    for i in range(il+1,ir-1):
        dn[i] = (n[i-1] - n[i+1])/(2.0*dx)
        dd[i] = (D[i-1] - D[i+1])/(2.0*dx)
    for i in range(il+2,ir-2):
        d2n[i] = (n[i-2] + n[i+2] - 2.0*n[i])/(4.0*dx**2)
    for i in range(il,ir):
        term1[i] = dn[i]*dd[i] - (dn[i+1]*dd[i+1] + dn[i-1]*dd[i-1])*0.5
        term2[i] = D[i]*d2n[i] - 0.5*(D[i-1]*d2n[i-1] + D[i+1]*d2n[i+1])*0.5
    #term1 = dn*dd
    #term2 = D*d2n
    dn_out = -(term1 + term2)*dt
    n_out = n + dn_out
    flux = 1.0
    return flux,n_out
    
def diff_F2(r,n,t, r_cent,ind_low,D):
    
    #get diffusion properties    
    dr = r[1] - r[0]
    Z = (r - r_cent)/(np.sqrt(t*D*4.0))
    f = (np.exp(-Z**2))/(np.sqrt(4.0*np.pi*D*t))
    n_diff = scisig.convolve(n, f, mode = 'same')
    n_diff*= dr
    #j = ro.romberg_samp(n_diff,r)
    #print('Mass after diffusion = ' + str(j))
    #ax.plot(x,B*n_diff/j, label = 't = ' + str(c))
    #check = ro.romberg_samp(B*n_diff/j,r)
    #print('Check is ' + str(check) + ' at time ' + str(t))
    return n_diff#/j