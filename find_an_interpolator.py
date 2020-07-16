# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:32:16 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 

import scipy.interpolate as spint 

import os 

direc = os.getcwd() 

file = np.loadtxt(os.path.join(direc, 'many_iteration','poisson','1500eV',
                               'analysed_outputs','19','density_t14559.txt'))

r_file = file[:,1]
d_file = file[:,0]

r_file = np.flip(r_file, axis = 0)
d_file = np.flip(d_file, axis = 0)

extra = 0.001
r = np.linspace(r_file[0]+ extra, r_file[-1], num = 1000)

d_pchip = spint.pchip_interpolate(r_file, d_file,r)
d_krogh = spint.krogh_interpolate(r_file, d_file, r)
d_bary = spint.barycentric_interpolate(r_file,d_file,r)
fig, ax = plt.subplots()
ax.plot(r_file, d_file, label = 'original')
ax.plot(r, d_pchip, label = 'pchip')
ax.plot(r,d_krogh, label = 'krogh')
ax.plot(r,d_bary, label = 'barycentric')
plt.yscale('log')
plt.legend()
plt.show()
