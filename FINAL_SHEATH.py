# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:16:51 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import rkf45_sheath_FINAL as rkf 
import scipy.interpolate as spint

"""Real units"""
e = 1.6*10**(-19) # in C
kb = 1.38*10**(-23) # in J/K
mi = 1.67*10**(-27) #in kg
n0 = 10**19 #in m^-3
eps0 = 8.854*10**(-12) # in F/m
me = 9.11*10**(-31) # in kg
phi_0 = 1000.0 # in V


T = phi_0*e/kb # in K
u0 = np.sqrt(kb*T/mi) # in m/s
ld = (eps0*kb*T/(e*e*n0))**(0.5) # in m
epsp = mi*u0*u0/(2.0*e) #in V


"""Normalised Units"""

T_n = T*kb/e #in V
T_n /= phi_0 #non-dim
epsp_n = 1.01*epsp/phi_0 #non-dim
e_n = e/e #non-dim
n_n = n0/n0 # non-dim

#FInal constant
K = (2.0*n0*e*ld*ld/(phi_0*eps0))


h = 0.01
xl = 0.0
xr = 100.0

phi,x, h_arr, err_arr = [],[],[],[]
phi_f = -0.0001
phi_s = np.log((2.0*mi/(2.0*np.pi*me))**0.5)
phi = [-phi_s - 0.5] # factor of 0.5 comes from pre-sheath addition

x = [0.0]
def sheath(x,phi_n):
    term1 = T_n*np.exp(phi_n/T_n) - T_n
    term2 = 2.0*epsp_n*(1.0 - phi_n/epsp_n)**0.5 - 2.0*epsp_n
    diff = K*(term1 + term2)
    diff = diff**0.5
    return diff
count = 0
while phi[-1] < phi_f:
    if count ==137:
        print( count)
    x_out, phi_out, h, err = rkf.rkf(x[-1], phi[-1], h, sheath, xr,xl)
    x = np.append(x,x_out)
    phi = np.append(phi,phi_out)
    h_arr = np.append(h_arr, h)
    err_arr = np.append(err_arr, err)
    count +=1

fig, ax  = plt.subplots()
ax.plot(x,phi)
plt.show()
