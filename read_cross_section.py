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

cs_ion = spint.pchip_interpolate(ion[:,0],ion[:,1], en)
cs_eff = spint.pchip_interpolate(eff[:,0], eff[:,1], en)

r_rate = cs_ion/cs_eff

#np.savetxt(direc + os.sep + 'h2_ion_rate.txt', r_rate)
with open(direc + os.sep + 'h2_ion_rate.txt', 'w') as f:
    f.write('%.15f' % r_rate )