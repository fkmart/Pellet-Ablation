# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 09:31:51 2020

@author: Kyle
"""
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

index = 5

l = 2000
step = 1
fig, ax  = plt.subplots()
ax.set_xlabel(r'$\tilde{t}$', fontsize = 12, rotation = 0 )
ax.set_ylabel(r'$R$', fontsize = 12, rotation = 0)
ax.xaxis.set_label_coords(0.52,-0.03)
ax.yaxis.set_label_coords(-0.04,0.42)

t_plot = np.linspace(0.1,0.2, num = l)
plot_ratio = np.zeros(l)
for j in range(1, l, step):
    load_dir = os.path.join(direc, 'many_iteration', 'sheath', '1000eV', 'analysed_outputs') + os.sep
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    dens = file[5:,1]
    dens_s = np.asarray([np.abs(float(w)) for w in dens])

    
    load_dir = os.path.join(direc, 'many_iteration', 'ion', '1000eV', 'analysed_outputs') + os.sep
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    
    dens = file[5:,1]
    dens_i = np.asarray([np.abs(float(w)) for w in dens])

    dens_ratio = dens_i[index]/dens_s[index]
    plot_ratio[j] = dens_ratio
ax.plot(t_plot[1:],plot_ratio[1:])#, label = r'$\Delta \tilde{t} = $' +'{:3.2f}'.format((j-1)*delta_t))
plt.yticks([0.00,0.50,1.00,1.50,2.00, 2.50])
#plt.yscale('log')
#plt.xlim(r_s[0], 6.5)
#plt.xscale('log')
#plt.legend(ncol = 2)
plt.savefig('elec_dens_ratio_r1.png', format = 'png', dpi = 1400)
plt.show()