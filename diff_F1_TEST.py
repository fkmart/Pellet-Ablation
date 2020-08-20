# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 09:36:16 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt
import diffusion_F1 as diff
import os 
from gen_var import diff_ion_sca, r0

direc = os.getcwd()

load_dir = os.path.join(direc,'many_iteration','neutral','500eV','analysed_outputs') + os.sep 

file = np.genfromtxt(load_dir + 'outputs1.txt', delimiter = ',', dtype = 'str')

dens = file[5:,1]
dens = np.asarray([float(w) for w in dens])
r = file[5:,0]
r = np.asarray([float(w) for w in r])

ind1 = 25
ind2 = 726

#get diffusion properties 
neut_dens = 0.01*(1.0 + r[25]**2)/(1.0 + r**2)
D = neut_dens**-1
#D = 50.0
t = 0.0005 

flux ,n_diff = diff.diff_flux(dens,r,D,ind1,ind2,t)
n_diff *= t*diff_ion_sca/r0

fig, ax = plt.subplots()
plt.yscale('log')
ax.plot(r,dens)
ax.plot(r,n_diff)
plt.show()