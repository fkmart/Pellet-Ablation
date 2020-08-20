# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 14:15:02 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 

direc = os.getcwd() 
lab = ['0.1','0.3','0.5','0.7','0.9']
c = 0
fig, ax = plt.subplots()
for t in range(50,500, 100):
    load_dir = os.path.join(direc, 'one_iteration','analysed_outputs' , 't_' + str(t)) + os.sep
    file = np.loadtxt(load_dir + 'norm_density.txt')
    tlab = t/500
    ax.plot(file[0,4:], np.flip(file[1,:-4],axis = 0), label = r'$\tilde{t} = $' + lab[c] )
    c+=1
plt.xlabel(r'$\tilde{r}$')
plt.ylabel(r'$\tilde{n}_e$', rotation = 0)
ax.xaxis.set_label_coords(0.5,-0.04)
ax.yaxis.set_label_coords(-0.04, 0.51)
plt.yscale('log')
plt.legend()
#plt.savefig('elec_density.png', format = 'png', dpi = 1400)
plt.show()