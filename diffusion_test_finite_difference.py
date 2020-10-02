# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 13:42:36 2020

@author: Kyle
"""
import numpy as np
import matplotlib.pyplot as plt
from gen_var import delta_t, cs_ion_neut, v_ion, rp_hr, pel_dens_numb, r0, tf
import os 

def solver2(n_old, D, h, delta_t, n_new):
    for i in range(1, len(n_old) - 1):
        term1 = D[i] *(n_old[i+1] + n_old[i-1] - 2.0*n_old[i])/(h*h)
        term2 = ((D[i+1] - D[i])/h)*(n_old[i+1] - n_old[i])/h
        n_new[i] = n_old[i] + delta_t *(term1 + term2)
    return n_new

direc = os.getcwd()

load_dir = os.path.join(direc,'many_iteration','neutral', '1000eV','analysed_outputs') + os.sep
file = np.genfromtxt(load_dir +  'outputs1.txt', dtype = 'str', delimiter = ',')

e_dens = file[5:,1]
r = file[5:,0]

r_arr = np.zeros(len(r))
dens_arr = np.zeros(len(e_dens))
for a in range(0 , len(r)):
    dens_arr[a] = float(e_dens[a])
    r_arr[a] = float(r[a])
e_dens = dens_arr[4:94]
r = np.asarray(r_arr[4:94]) 

t_ind = 100000
dens = pel_dens_numb*0.01*(1 + rp_hr[t_ind])/(1 + r**2)
mfp = 1.0/(dens*cs_ion_neut)
D = mfp *v_ion
D_norm = D*tf/(r0*r0)

h = r[1] - r[0]

n_out = np.zeros(len(r))
fig, ax = plt.subplots()
ax.plot(r,e_dens, label = 'initial')
coeff = D_norm*(delta_t/10.0)/(h**2)
print('maximum coefficient is ' + str(np.max(coeff)))
for a in range(0, 10):
    n_out = solver2(e_dens, D_norm, h, delta_t/10.0,n_out)
    if (a%2 == 1):
        ax.plot(r,n_out,label = r'$\Delta t = $' + str(a+1) )
    e_dens = np.copy(n_out)
plt.yscale('log')
plt.legend()
plt.show