# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 14:12:24 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 

direc = os.getcwd() 
TYPE = ['neutral','poisson','sheath']

length = 20
t = np.linspace(1,19,num = 19)
fig,ax = plt.subplots()
n_loss = np.zeros(length -1)
for i in range(0, len(TYPE)):
    for j in range(1,length):
        load_dir = os.path.join(direc,'many_iteration',TYPE[i],'1000eV','analysed_outputs') + os.sep
        file = np.genfromtxt(load_dir + 'outputs' + str(j)+'.txt', dtype = 'str', delimiter = ',')
        n_loss[j-1] = float(file[3,4])
    ax.plot(t,n_loss)
plt.legend()
plt.yscale('log')
plt.show()