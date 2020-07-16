#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 12:45:54 2019

@author: kyle
"""
#NEW CHANGES HAVE BEEN MADE NOW THAT WHILE LOOP ACTS SUCCESSFULLY, EACH ENERGY HAS A DIFFERENT NUMBER OF ELEMENTS IN OUTPUT
#Module to determine particle energies at pellet surface

import numpy as np 
import stop_analysis_stop_point
import os
def term_energy(particle,t, le, direc):
    from electron import RME, le
    from gen_var import I, rp, style 
    if style == 'many':
        from gen_var import rp_hr as rp
    else:
        pass
    arrsave1 = []
    ind= []
    for j in range(0,le -1):
        
        ener_arr = np.load(os.path.join(direc,'EvsR_' + particle +'_E0'+str(j)+'.npy'))
       
        """
        ind.append(len(ener_arr))

        if (np.ma.size(ener_arr, axis = 0)<lr-1):
            pass
        else:       
            
            arrsave1.append(ener_arr[0,0])
            arrsave1.append(ener_arr[lr-2,0])
        #arrsave2 = stop_analysis_stop_point.stop_point(ener_arr, ind, arrsave2, j, r)"""

        if (ener_arr[-1,1] < rp[t]) or ener_arr[-1,1] == rp[t]:
            arrsave1.append(ener_arr[0,0])
            arrsave1.append(ener_arr[-1,0])
        else:
            pass
    arrsave1 = np.reshape(arrsave1, (int((len(arrsave1))/2),2), order = 'C')
    
    return arrsave1, ind#, arrsave2

