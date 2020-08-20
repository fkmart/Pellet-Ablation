# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:30:34 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import many_start

types = ['neutral','poisson', 'sheath', 'ion']

k = 4

fig, ax = plt.subplots()
direc = os.getcwd()
for i in range(0,len(types)):
    load_dir = os.path.join(direc,'many_iteration', types[i], 'analysed_data') + os.sep
    shifts = np.loadtxt(load_dir + 'shifts.txt')
    time = many_start + int(shifts[k])
    file = np.loadtxt(load_dir + 't_' + time + os.sep + 'accummulated_density.txt')
    ax.plot(file[:,1], file[:,0], label = types[i])
plt.show()
        