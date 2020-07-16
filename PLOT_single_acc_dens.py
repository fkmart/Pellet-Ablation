# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 14:57:17 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 

import os 

direc = os.getcwd() 
fig, ax = plt.subplots()
TYPE = ['neutral','poisson','sheath']
ls = ['-','--',':']
ax.set_xlim(0.0,3)
ax.set_ylim(0.00001,0.0001)
ax2 = ax.twinx()

for k in range(0,len(TYPE)):
     load_dir = os.path.join(direc, 'many_iteration', TYPE[k], '1500eV','analysed_outputs', '0')+os.sep
     file = np.loadtxt(load_dir + 'accummulated_density_t12792.txt')
     file2 = np.loadtxt(load_dir + 'potential_t12792.txt')
     ax.plot(file[1,:], file[0,:], label = TYPE[k], color = 'green', linestyle = ls[k])
     ax2.plot(file2[1,:], file2[0,:], label = TYPE[k], color = 'blue', linestyle = ls[k])
plt.legend()
ax.set_yscale('log')
plt.show()