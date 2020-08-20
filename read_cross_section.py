# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 09:59:52 2020

@author: Kyle
"""

import os
import numpy as np
import scipy.interpolate as spint 
from gen_var import trunc_fac, I

direc = os.getcwd() 

en = trunc_fac*I
ion = np.loadtxt(direc + os.sep +'h2_ionisation.txt')
eff = np.loadtxt(direc + os.sep + 'h2_effective.txt')
elas = np.loadtxt(direc + os.sep + 'h2_elastic.txt')

cs_ion = spint.pchip_interpolate(ion[:,0],ion[:,1], en)
cs_eff = spint.pchip_interpolate(eff[:,0], eff[:,1], en)
cs_elas = spint.pchip_interpolate(elas[:,0], elas[:,1], en)

r_rate = cs_ion/cs_eff

#np.savetxt(direc + os.sep + 'h2_ion_rate.txt', r_rate)
with open(direc + os.sep + 'h2_ion_rate.txt', 'w') as f:
    f.write('%.15f' % r_rate )
    
#now check that remaining probability of interaction comes from excitation
excitation_cs = np.zeros(15)
for i in range(1,16):
    file = np.loadtxt(direc + os.sep + 'h2_excitation' + str(i) + '.txt')
    excitation_cs[i-1] = spint.pchip_interpolate(file[:,0], file[:,1], en)

cs_exc = np.sum(excitation_cs)

cs_total = cs_exc + cs_ion #+ cs_elas
