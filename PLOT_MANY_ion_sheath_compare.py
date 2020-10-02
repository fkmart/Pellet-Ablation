# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 18:51:21 2020

@author: Kyle
"""
import os 
import matplotlib.pyplot as plt 
import numpy as np 
from gen_var import delta_t

direc = os.getcwd() 
TYPE = 'ion'

l = 2000
step = 400
fig, ax  = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0 )
ax.set_ylabel(r'$R$', fontsize = 12, rotation = 0)
ax.xaxis.set_label_coords(0.45,-0.03)
ax.yaxis.set_label_coords(-0.04,0.46)

for j in range(1, l, step):
    load_dir = os.path.join(direc, 'many_iteration', 'sheath', '1000eV', 'analysed_outputs') + os.sep
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    dens = file[5:,1]
    dens_s = np.asarray([np.abs(float(w)) for w in dens])
    r = file[5:,0]
    r_s = np.asarray([float(w) for w in r])
    
    load_dir = os.path.join(direc, 'many_iteration', 'ion', '1000eV', 'analysed_outputs') + os.sep
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    
    dens = file[5:,1]
    dens_i = np.asarray([np.abs(float(w)) for w in dens])
    r = file[5:,0]
    r_i = np.asarray([float(w) for w in r])
    
    starter = 40
    low_s = next(p[0] for p in enumerate(dens_s) if p[1] > 0.0)
    up_s = next(p[0] for p in enumerate(dens_s[starter:]) if p[1] == 0.0)
    
    low_i = next(p[0] for p in enumerate(dens_i) if p[1] > 0.0)
    up_i = next(p[0] for p in enumerate(dens_i[starter:]) if p[1] == 0.0)
    
    print(low_s - low_i)
    print(up_s - up_i)
    if low_s < low_i:
        low = low_s
    else:
        low = low_i
    
    if up_i < up_s:
        up = up_i
    else:
        up = up_s
    up +=starter 
    low +=1
    dens_ratio = dens_i[low:up]/dens_s[low:up]
    ax.plot(r_s[low:up],dens_ratio, label = r'$\Delta \tilde{t} = $' +'{:3.2f}'.format((j-1)*delta_t))
#plt.yscale('log')
plt.xlim(r_s[0], 6.5)
#plt.xscale('log')
plt.legend(ncol = 2)
plt.savefig('elec_dens_ratio_all_t.png', format = 'png', dpi = 1400)
plt.show()