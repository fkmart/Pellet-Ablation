# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 10:07:50 2020

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
ax.set_ylabel(r'$\tilde{n}_i$', fontsize = 12, rotation = 0)
ax.xaxis.set_label_coords(0.45,-0.03)
ax.yaxis.set_label_coords(-0.04,0.46)

for j in range(1, l, step):
    load_dir = os.path.join(direc, 'many_iteration', TYPE, '1000eV', 'analysed_outputs') + os.sep
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    dens = file[5:,2]
    dens = np.asarray([float(w) for w in dens])
    r = file[5:,0]
    r = np.asarray([float(w) for w in r])
    ax.plot(r,dens, label = r'$\Delta \tilde{t} = $' +'{:3.2f}'.format((j-1)*delta_t))
plt.yscale('log')
plt.xlim(r[0], 6.5)
#plt.xscale('log')
plt.legend(ncol = 2)
plt.savefig('acc_ion_dens_' + TYPE +'_log.png', format = 'png', dpi = 1400)
plt.show()