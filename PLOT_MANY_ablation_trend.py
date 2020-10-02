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

TYPE = ['neutral', 'poisson', 'sheath', 'ion']

limit = 2000
jump = 1
rpd_plot = np.zeros(int(limit/jump))
t_plot = np.arange(0,limit,jump)
t_plot = np.linspace(many_start/10**6,(many_start + 50*2000)/10**6, limit)

rpd_evap = rpd_hr[many_start:many_start + 2000*50:50*jump]
rpd_NGS = t_plot**(-2.0/5.0)
fig, ax = plt.subplots()
ax.plot(t_plot, np.abs(rpd_evap)/np.amax(rpd_evap), label = 'Evaporative')
ax.plot(t_plot, np.abs(rpd_NGS)/np.amax(rpd_NGS), label = 'NGS')
plt.yscale('log')
plt.ylabel(r'$\dot{\tilde{r}}_p$', fontsize = 12, rotation = 0)
plt.xlabel(r'$\tilde{t}$')
ax.xaxis.set_label_coords(0.55,-0.06)
ax.yaxis.set_label_coords(-0.04,0.42)
labels = ['Neutral','Charged Cloud','Sheath', 'Diffusive Ion']
for t in range(0, len(TYPE)):
    load_dir = os.path.join(direc, 'many_iteration',TYPE[t],'1000eV','analysed_outputs') + os.sep
    for j in range(1, limit,jump):
        file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
        N_loss = float(file[3,4])
        current_i = int(file[1,4])
        new_rp = (rp_hr[current_i]**3 - N_loss)**(1.0/3.0)
        drp = rp_hr[current_i] - new_rp
        rpd = drp/(delta_t)
        rpd_plot[j] = rpd
    ax.plot(t_plot,np.abs(rpd_plot)/(np.amax(rpd_plot)), label = labels[t])
plt.legend()
plt.savefig('rpd_trends.png', format = 'png', dpi = 1400)
plt.show()