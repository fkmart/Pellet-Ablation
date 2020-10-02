# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 10:36:57 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint 
from gen_var import balmer_rate, rp_hr, rc_hr, eps, pel_dens_numb, dens_plas, delta_t
import os 

direc = os.getcwd() 

TYPE = 'neutral'

l = 2000
step = 400
fig, ax  = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0 )
ax.set_ylabel('Intensity/ Arbitrary units', fontsize = 10, rotation = 90)
ax.xaxis.set_label_coords(0.45,-0.03)
ax.yaxis.set_label_coords(-0.10,0.52)

for j in range(1, l, step):
    load_dir = os.path.join(direc, 'many_iteration', TYPE, '1000eV', 'analysed_outputs') + os.sep
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
    ax.plot(r,emission, label = r'$\Delta \tilde{t} = $' +'{:3.2f}'.format((j-1)*delta_t))
plt.yscale('log')
plt.xlim(r[0], 8.5)
#plt.xscale('log')
plt.legend(ncol = 3, loc = 'lower center')
#plt.savefig('potential_' + TYPE +'_log.png', format = 'png', dpi = 1400)
plt.show()