# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:50:59 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 
import iterative_sol as SOR 
from gen_var import rp,rc , r0, dens_plas, e, epsilon0
import discret_mat as dm 
import scipy.interpolate as spint 
from electron import RME,M_fac
import norming_density as nd


direc = os.getcwd()
load_dir = os.path.join(direc, 'many_iteration', 'poisson','analysed_outputs','t_20') + os.sep
file = np.loadtxt(load_dir + 'density_t20.txt')

pos = file[:,1]
bins = file[:,0]

g = spint.interp1d(pos,bins, kind = 'cubic', fill_value = 'extrapolate')

r = np.linspace(rp[20], rc[20], num = 513, endpoint = 'true')
phi = np.zeros(513)
phi[0] = -2800/(RME*M_fac)
non_dim = r0*r0*e*dens_plas/(epsilon0*RME*M_fac)
A = dm.discret(r)
bins_interp = g(r)
#bins_interp = spint.pchip(pos, bins, r, extrapolate = 'true')

e_dens = bins_interp/(np.sum(bins_interp))
phi_out = SOR.SOR(A,phi, non_dim*e_dens, r)
fig, ax = plt.subplots()
ax.plot(pos,bins)
ax.plot(r,bins_interp, color = 'purple', linestyle = '--')
ax.plot(r,e_dens, color = 'red', linestyle = '--')
plt.yscale('log')
ax2= ax.twinx()
ax2.plot(r,M_fac*RME*phi_out, color = 'orange')
plt.show()