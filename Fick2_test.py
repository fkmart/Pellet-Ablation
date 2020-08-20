# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:42:36 2020

@author: Kyle
"""

import numpy as np
import os 
import matplotlib.pyplot as plt
from diffusion_F1 import diff_F2 as f2
import romberg as ro
from gen_var import h2_ion_rate, eps, rp_hr, delta_t, ratio, ion_dens_renorm
from gen_var import v_ion, pel_dens_numb, cs_ion_neut, r0, tf

direc = os.getcwd()

load_dir = os.path.join(direc,'many_iteration','neutral','500eV','analysed_outputs') + os.sep 

file = np.genfromtxt(load_dir + 'outputs1.txt', delimiter = ',', dtype = 'str')

dens = file[5:,1]
dens = np.asarray([float(w) for w in dens])
r = file[5:,0]
r = np.asarray([float(w) for w in r])

ind1 = 26
ind2 = 759

i = 65792
ind_max = np.argmax(dens)
r_cent = r[ind_max]
r_cent = 0.5*(r[-1] + r[0])
#get diffusion properties 
neut_dens = 0.01*(1.0 + r[25]**2)/(1.0 + r**2)
dens_neut_n = np.zeros(len(dens))
diff_ion_n = np.zeros(len(dens))
D = neut_dens**-1
t = delta_t/ratio

n_ion = dens*h2_ion_rate
dens += n_ion # updates electron density with the released ionised electrons
dens_neut_n[ind1:ind2] = eps*(1.0 + rp_hr[i]**2)/(1.0 + r[ind1:ind2]**2)
diff_ion_n[ind1:ind2] = dens_neut_n[ind1:ind2]**-1 
diff_ion_n[diff_ion_n == 0.0] = diff_ion_n[ind1]
r_cent = 0.5*(r[-1] + r[0])

#scaling quantities
D_sca = v_ion/(pel_dens_numb*cs_ion_neut)
n_new_sca = (np.exp(r0*r0/(D_sca*tf)))*(1.0/(np.sqrt(D_sca*tf)))


fig,ax = plt.subplots()
ax.plot(r,n_ion)
count = 1
dens[ind1] = 0.0
#B = ro.romberg_samp(dens,r)
ablating_ions = 0.0
print(np.sum(n_ion))
while count < 5:
    B = np.sum(n_ion)
    n_ion = f2(r,n_ion,t, r_cent, ind1-1,D)
    C = np.sum(n_ion)
    n_ion *= B/np.sum(n_ion)
    print(np.sum(n_ion))
    n_ion[ind2:] = 0.0
    ablating_ions += np.sum(n_ion[:ind1])
    ax.plot(r,n_ion)
    n_ion[:ind1] = 0.0
    #t += 0.000005
    count += 1
#ablating_ions *= ion_dens_renorm
plt.legend()
plt.yscale('log')
plt.show()