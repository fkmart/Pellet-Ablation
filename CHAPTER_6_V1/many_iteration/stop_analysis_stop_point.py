#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:35:52 2019

@author: kyle
"""

import numpy as np
from gen_var import *
import os
def stop_point(arr,ind,particle, r, i, le_mid, direc):
    lr = len(r)
    arrsave2 = []
    for j in range(0,le_mid):
        ener_arr = np.load(os.path.join(direc,'EvsR_' + particle + '_t'+str(i) +'_E0'+str(j)+'.npy'))
        arrsave2.append(ener_arr[0,0])
        x = ind[j]
        arrsave2.append(r[lr -1 - x])
    arrsave2 = np.reshape(arrsave2, (int(len(arrsave2)/2),2), 'C')
    #np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_deposition_cloud_1keV_fs_t' + str(i) +'.txt'), arrsave2, fmt = ('%f'))       
    return arrsave2