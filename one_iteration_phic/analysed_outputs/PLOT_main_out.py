#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 10:44:08 2019

@author: kyle
"""
#plotting file to check new outputs

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spint
#from gen_var import dt, pel_pot, lp 
import os
import scipy.integrate as spinteg

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

particle = 'electron'
print('Particle type is ' + str(particle))


i = 50

for p in range(200,2000, 200):
    file = np.loadtxt('real_density_pot_test_t'+str(i) +'pot-'+str(float(p))+'.txt')
    #file = np.loadtxt(GNU_'+particle+'_deposition_cloud_1keV_fs_t' + str(i) + '.txt')
    
    #file = np.loadtxt('GNU_'+particle+'_terminal_ener_1keV_fs_t' + str(i) + '.txt')
    ax.plot(file[:,0], file[:,1], label = r'$\phi = '+str(p)+'$')


#file = np.loadtxt(os.path.join('Analysed Outputs', 'GNU_'+particle+'_min_imp_ener_1keV_fs_t.txt'))

#t1 = file[:,1]
#t_int =  np.arange(0.0,0.999,0.001)
#g = spint.interp1d(t1,file[:,0], kind = 'cubic' )
#g= spint.pchip_interpolate(t1, file[:,0], t_int)
#t2 = np.linspace(t1[0],t1[-1], 90)
#term_en = g(t2)
#ax.plot(t1, file[:,0], 'x')#, label = r'$\frac{t}{t_f} = 0.'+str(i)+ '$')
#ax.legend(ncol = 2, fontsize = 8.5)#, loc = 3)
ax.set_ylabel(r'$\frac{\rho_e}{\rho_t}$', fontsize = 15, rotation = 0)
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.xaxis.set_label_coords(0.53,-0.03)
ax.yaxis.set_label_coords(-0.08, 0.440)
ax.xaxis.labelpad = -10
ax.yaxis.labelpad = -2.5
#ax.set_title('Minimum impacting energy of electrons on pellet surface with time')
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.legend()
plt.show()
#ax.set_xscale('log')
#plt.savefig('stopped_electron_dens.png', dpi = 800, format = 'png')

    
    