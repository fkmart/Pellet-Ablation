# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:19:06 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp_hr, eps, pel_dens_numb, dens_plas, balmer_rate, delta_t

direc = os.getcwd() 

fig, ax  = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0 )
ax.set_ylabel('Intensity/ Arbitrary units', fontsize = 10, rotation = 90)
ax.xaxis.set_label_coords(0.50,-0.03)
ax.yaxis.set_label_coords(-0.10,0.52)

TYPE = ['neutral','poisson','sheath','ion']

for i in range(0, len(TYPE)):
    load_dir = os.path.join(direc, 'many_iteration', TYPE[i], '1000eV', 'analysed_outputs') + os.sep
    j = 800
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    e_dens = file[5:,1]
    e_dens = np.asarray([float(w) for w in e_dens])
    t_index = int(file[1,4])
    r = file[5:,0]
    r = np.asarray([float(w) for w in r])
    neut_dens = eps*(1.0 + rp_hr[t_index]**2)/(1.0 + r**2)
    e_dens *= dens_plas/pel_dens_numb
    emission = balmer_rate*(e_dens*neut_dens)/(e_dens + neut_dens)
    emission /= np.amax(emission)
    ax.plot(r,emission, label = TYPE[i])
plt.yscale('log')
plt.xlim(r[0], 8.5)
plt.legend(ncol = 2, loc = 'lower left')
plt.savefig('balmer_emission_comparison.png', format = 'png', dpi = 1400)
plt.show()