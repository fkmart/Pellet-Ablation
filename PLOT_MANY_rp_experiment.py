# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:55:18 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp_hr, tf, delta_t, r0, N_0_sca, rpd_hr, many_start

direc= os.getcwd() 

TYPE = ['poisson', 'sheath', 'ion']

limit = 2001
jump = 1
rp_plot = np.zeros(int(limit/jump))
rp_plot[0] = rp_hr[many_start]
t_plot = np.arange(0,limit,jump)
t_plot = np.linspace(many_start/10**6,(many_start + 50*(limit-1))/10**6, limit)


fig, ax = plt.subplots()
plt.yscale('log')
plt.ylabel(r'$\tilde{r}_p$', fontsize = 12, rotation = 0)
plt.xlabel(r'$\tilde{t}$', fontsize = 12)
ax.xaxis.set_label_coords(0.55,-0.06)
ax.yaxis.set_label_coords(-0.04,0.42)
labels = ['Charged Cloud','Sheath', 'Diffusive Ion']
for t in range(0, len(TYPE)):
    load_dir = os.path.join(direc, 'many_iteration',TYPE[t],'1000eV','analysed_outputs') + os.sep
    for j in range(1, limit,jump):
        file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
        current_i = int(file[1,4])
        rp_plot[j] = rp_hr[current_i]
    ax.plot(t_plot,rp_plot, label = labels[t])
plt.legend()
plt.savefig('rp_experiment.png', format = 'png', dpi = 1400)
plt.show()