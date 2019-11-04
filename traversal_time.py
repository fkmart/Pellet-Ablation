#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 15:36:15 2019

@author: kyle
"""

import numpy as np
import scipy.integrate as spit


particle  = 'electron'

me = 9.11*10**(-31)
e = 1.6*10**(-19)
t_final = []
for i in range(1,10):
    print(i)
    t= []
    
    for j in range(0, 398):
        en = np.load('EvsR_'+particle +'_t%d_E0%d.npy' % (i,j))
        integrand = np.sqrt(2.0*e*en[:,0]/me)
        x = en[:,1]**10**(-3)
        t1 = spit.simps(integrand**(-1), en[:,1])
        t.append(t1)
    t_final.append(t)

np.reshape(t_final, (398,9,), 'F')
        
        