# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:33:39 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 
import scipy.interpolate as spint 
from gen_var import rp, rc


direc = os.getcwd() 
load_dir = os.path.join(direc, 'many_iteration', 'poisson','analysed_outputs','t_20') + os.sep
file = np.loadtxt(load_dir + 'density_t20.txt')

r1 = np.linspace(rp[20], rc[20], num = 512, endpoint = 'true')
r2 = np.linspace(rp[21], rc[21], num = 512, endpoint = 'true')

g = spint.interp1d(file[:,1], file[:,0], kind = 'cubic', fill_value = 'extrapolate')

f = g(r1)

fig, ax  = plt.subplots() 
ax.plot(file[:,1], file[:,0])
ax.plot(r1,f)
plt.yscale('log')
plt.show()