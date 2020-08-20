# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 12:30:21 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 

direc = os.getcwd() 
TYPES = ['neutral','poisson','sheath','ion']
energies = ['500eV','1000eV']

c = ['blue','red','green']

ls = ['-','--','-.',':']

for k in range(1,200, 20):
    fig,ax = plt.subplots(3,1)
    for i in range(0,len(TYPES)):
        for j in range(0, len(energies)):
            file = np.genfromtxt(os.path.join(direc,'many_iteration',TYPES[i], energies[j], 'analysed_outputs', 'outputs' + str(k) + '.txt'), delimiter = ',', dtype = 'str')
            r = file[5:,0]
            r = np.asarray(r, dtype = 'float')
            e_dens = file[5:,1]
            e_dens = np.asarray(e_dens, dtype = 'float')
            pot = file[5:,3]
            pot = np.asarray(pot, dtype = 'float')
            ax[j].plot(r,e_dens, linestyle = ls[i], color = c[j])
    plt.show()