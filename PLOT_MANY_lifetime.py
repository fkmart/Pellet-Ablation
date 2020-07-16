# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 09:23:08 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 
from gen_var import rgl, ener_ablat_sca
import os

direc = os.getcwd()
TYPE = ['neutral','poisson','sheath']
style = ['-','--',':']
col = ['blue','red','green']
energy = ['500eV','1000eV','1500eV']
fig,ax = plt.subplots()
#plt.ylim(10**13,10**18)
t_plot = np.linspace(0,19, num = 20)
for j in range(0,len(TYPE)):
    for k in range(0,len(energy)):
        load_dir = os.path.join(direc, 'many_iteration',TYPE[j], energy[k],'analysed_outputs') + os.sep
        indices = np.loadtxt(load_dir + 'time_indices.txt')
        #t_params = np.loadtxt(load_dir + 't_params.txt')
        imp_frac = np.loadtxt(load_dir + 'impacting_fractions.txt')
        lifetimes = np.loadtxt(load_dir + 'pellet_lifetime.txt')
        ret_flux = np.loadtxt(load_dir + 'retarded_ener_flux.txt')
        ax.plot(t_plot,lifetimes[:,0], color = col[k], linestyle = style[j], label = energy[k]+'-'+TYPE[j] )
plt.yscale('log')
#plt.legend(ncol = 3)
plt.xlabel(r'$\Delta t/ 5\mu \mathrm{s}$')
plt.ylabel(r'$\tilde{t}_f/ \Delta t$', rotation = 0)
plt.savefig('lifetimes.png', format = 'png', dpi = 1200)
plt.show()