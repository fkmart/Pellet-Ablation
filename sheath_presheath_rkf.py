# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 10:12:25 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import rkf45_sheath_presheath as rkf 
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
K = (2.0*n0*e*ld*ld/(phi_0*eps0))**(-0.5)

def sheath_func(phi_n,x):
    term1 = T_n*np.exp(phi_n/T_n) - T_n*np.exp(0.5)
    term2 = 2.0*epsp_n*(1.0 - phi_n/epsp_n)**0.5
    term3 = -2.0*epsp_n*(1.0 - T_n/(2.0*epsp_n)**0.5)
    diff1 = term1 + term2 + term3
    diff = diff1**(-0.5)
    diff*= K**-1
    return diff

phi_s = np.log(np.sqrt(2.0*mi/(2.0*np.pi*me)))
print(phi_s)
#l = 3000
phi_1 = 0.5 # in normalised
fac = 0.999
#phi = -np.linspace(phi_s, phi_1, l)
#phi = np.logspace( np.log10(phi_s),np.log10(phi_1), l)
x0 = 0.0
h = 0.01

#epsp_n *=phi_s*1.01

phi, x, dr_arr, err_arr = [],[],[],[]
phi = [-phi_s]
x = [x0]
count =0
while phi[-1] < phi_1*fac:
    if count == 314:
        print('here')
    phi_out, x_out, dr_out, err_out = rkf.rkf(phi[-1],x[-1],h,sheath_func, phi_1,-phi_s)
    phi = np.append(phi, phi_out)
    x = np.append(x, x_out)
    dr_arr = np.append(dr_arr, dr_out)
    err_arr = np.append(err_arr, err_out)
    count +=1
fig, ax= plt.subplots()
ax.plot(x, phi, color = 'violet')
#plt.yscale('log')
plt.show()







