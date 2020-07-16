# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 15:44:01 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 
import scipy.interpolate as spint

direc = os.getcwd() 

load_dir = os.path.join(direc, 'many_iteration', 'neutral','analysed_outputs','t_8780') + os.sep
file = np.loadtxt(load_dir + 'accummulated_density_t8780.txt')

pos = file[1,:]
bins = file[0,:]

indices = []
for i in range(1, len(bins) - 5):
    if bins[i] > 1.05*bins[i+1]:
        indices = np.append(indices,i)

bins_new = np.delete(bins,indices)
pos_new = np.delete(pos, indices)

f = spint.interp1d(pos_new, bins_new, kind = 'cubic')
bins_new = f(pos)
fig,ax = plt.subplots()
ax.plot(pos,bins)
ax.plot(pos, bins_new)
plt.show()