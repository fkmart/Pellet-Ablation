# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 11:47:33 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 

direc = os.getcwd() 
TYPE = ['neutral','poisson','sheath']

length = 2000
t = np.linspace(1, length -1,num = length - 1)
fig,ax = plt.subplots()
ax2 = ax.twinx()
ax2.plot(t,t**(-2.0/3.0), label = 'NGS Model', color = 'k')
ax2.set_ylabel('Arbitrary Units', fontsize = 10)
n_loss = np.zeros(length -1)
for i in range(0, len(TYPE)):
    for j in range(1,length):
        load_dir = os.path.join(direc,'many_iteration',TYPE[i],'1000eV','analysed_outputs') + os.sep
        file = np.genfromtxt(load_dir + 'outputs' + str(j)+'.txt', dtype = 'str', delimiter = ',')
        n_loss[j-1] = float(file[1,5])
    ax.plot(t,n_loss, label = TYPE[i])
plt.legend()
#plt.yscale('log')
plt.show()