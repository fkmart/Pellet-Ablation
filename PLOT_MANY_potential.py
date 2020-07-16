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
for i in range(0,20):
    fig,ax = plt.subplots(3)
    for k in range(0,len(energy)):
        for j in range(0,len(TYPE)):
            load_dir = os.path.join(direc, 'many_iteration',TYPE[j], energy[k],'analysed_outputs') + os.sep
            indices = np.loadtxt(load_dir + 'time_indices.txt')
            #t_params = np.loadtxt(load_dir + 't_params.txt')
            imp_frac = np.loadtxt(load_dir + 'impacting_fractions.txt')
            lifetimes = np.loadtxt(load_dir + 'pellet_lifetime.txt')
            ret_flux = np.loadtxt(load_dir + 'retarded_ener_flux.txt')
            file = np.loadtxt(load_dir + str(i) +os.sep+'potential_t'+str(int(indices[i])) +'.txt')
            l = ax[k].plot(file[1,:],np.abs(file[0,:]), color = col[k], linestyle = style[j], label = energy[k]+'-'+TYPE[j] )
            lab = energy[k]  
        ax[k].legend(l, [lab], loc = 'upper right')
        ax[k].set_yscale('log')
        ax[k].set_ylabel(r'$\left|\tilde{\phi}\right|$', rotation = 0)
        ax[k].set_xlabel(r'$\tilde{r}$')
        #ax[k].yaxis.set_label_coords(-0.02,0.0)
    plt.tight_layout()
    plt.savefig('many_potential_' + str(i) + '.png', dormat = 'png', dpi = 1400)
    plt.show()