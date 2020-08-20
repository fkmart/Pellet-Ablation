# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:00:28 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt
import os 
from gen_var import pel_pot

direc = os.getcwd() 

load_dir = os.path.join(direc, 'one_iteration_phic', 'analysed_outputs', 't_50') + os.sep

fig, ax = plt.subplots() 
for p in range(0, len(pel_pot),5):
    file = np.loadtxt(load_dir + 'lin_density_pot_' + str(pel_pot[p]) + '_test_mid.txt')
    #file = np.loadtxt(load_dir + 'lin_norm_dens_pot_' + str(pel_pot[p]) + '.txt')
    ax.plot(file[:,1], file[:,0], label = str(int(pel_pot[p])) + 'V')
    #ax.plot(file[0,:], file[1,:], label = str(int(pel_pot[p])) + 'V')
plt.legend()
plt.xlabel(r'$\tilde{r}$', rotation = 0, fontsize = 12)
plt.ylabel(r'$\tilde{n}_e$', rotation = 0, fontsize = 12)
ax.xaxis.set_label_coords(0.51,-0.04)
ax.yaxis.set_label_coords(-0.06, 0.54)
plt.savefig('faux_dens_lin_pot.png', format = 'png', dpi = 1400)
plt.show()