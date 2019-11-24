#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:35:52 2019

@author: kyle
"""

import numpy as np
from gen_var import *
import os
def stop_point(ind, r, p,a, le_mid, direc):
    lr = len(r)
    arrsave2 = []
    y = format(sig[a], '.2f')
    for j in range(0,le_mid):
        ener_arr = np.load(direc + 'prof' + str(j) +'_pot' + str(int(pel_pot[p])) + '_sig' +y+'.npy')
        arrsave2.append(ener_arr[0,0])
        x = ind[j]
        arrsave2.append(r[lr - x]) # removed a -1 
    arrsave2 = np.reshape(arrsave2, (int(len(arrsave2)/2),2), 'C')    
    return arrsave2