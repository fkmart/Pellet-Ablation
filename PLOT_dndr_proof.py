# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 11:29:36 2020

@author: Kyle
"""
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as spint
import os 

direc = os.getcwd() 

path = os.path.join(direc, 'one_iteration','proof', 'analysed_outputs','t_50') +  os.sep 
file = np.loadtxt(path + 'density.txt')

r1 = np.linspace(file[0,1], file[-1,1], num = 500, endpoint = 'true')
r2 = np.linspace(file[0,1], file[-1,1], num = 1000, endpoint = 'true')
r = np.flip(file[:,1], axis = 0)
d = np.flip(file[:,0], axis = 0)
d_tot = np.sum(d)
dens1 = spint.pchip_interpolate(r,d,r1)
dens1 *= d_tot/(np.sum(dens1))

dens2 = spint.pchip_interpolate(r,d,r2)
dens2 *= d_tot /(np.sum(dens2))
fig,ax = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0)
ax.set_ylabel(r'$\frac{\mathrm{d}n}{\mathrm{d}E}\Delta E \rightarrow G(r)$', fontsize = 8, rotation = 0)
ax.yaxis.set_label_coords(-0.08,0.58)
ax.scatter(file[:,1], file[:,0], label = 'raw data', marker = 'x', color = 'orange')
ax.plot(r1, dens1, label = 'grid1', color = 'navy')
ax.plot(r2, dens2, label = 'grid2', color = 'limegreen')
plt.yscale('log')
plt.legend()
plt.savefig('density_demo.png', format = 'png', dpi = 1200)
plt.show()