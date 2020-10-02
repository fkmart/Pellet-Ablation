# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:29:58 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 
from gen_var import many_start, rp_hr , t_hr, delta_t, N_0_sca

direc = os.getcwd() 
TYPE = 'sheath'

l = 2000
t = np.linspace(0,l,num = l)

fig, ax = plt.subplots()
ax2 = ax.twinx()

"NGS Baseline"

many_end = many_start + 2000*50
rpd_ngs = rp_hr[many_start:many_end]**(-2.0/3.0)
delta_r = rpd_ngs*delta_t 
proj_r = rp_hr[many_start:many_end] - delta_r
delta_N = rp_hr[many_start:many_end]**3 - proj_r**3 #Plot this and ignore scaling

dndt = np.zeros(l)

lines = ax2.plot(t, delta_N[::50], label = 'NGS Model', color = 'black')

load_dir = os.path.join(direc, 'many_iteration', TYPE, '1000eV', 'analysed_outputs') + os.sep
for j in range(1, l):
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', dtype = 'str', delimiter = ',')
    dndt[j] = file[3,4]
lines += ax.plot(t[:-1],dndt[1:], label = TYPE)

#DOING THE BEST FIT 

y_data = np.log10(dndt[1:])
fit = np.polyfit(np.log10(t[1:]), y_data, deg = 1)
y_fit = (10**(fit[1]))*t[1:]**fit[0]
lines += ax.plot(t[1:], y_fit, label = 'fit', color = 'red')


labs = [d.get_label() for d in lines] 
plt.yscale('log')
plt.legend(lines, labs)
plt.show()