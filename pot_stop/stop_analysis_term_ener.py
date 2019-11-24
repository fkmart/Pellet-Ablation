#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 12:45:54 2019

@author: kyle
"""
#NEW CHANGES HAVE BEEN MADE NOW THAT WHILE LOOP ACTS SUCCESSFULLY, EACH ENERGY HAS A DIFFERENT NUMBER OF ELEMENTS IN OUTPUT
#Module to determine particle energies at pellet surface

import numpy as np 
import os
def term_energy(r,p,a, le, direc):
    from electron import RME, le
    from gen_var import I, sig, pel_pot
    arrsave1 = []
    ind= []
    lr = len(r) 
    y = format(sig[a], '.2f')
    for j in range(0,le -1):
        ener_arr = np.load(direc + 'prof' + str(j) +'_pot' + str(int(pel_pot[p])) + '_sig' + y +'.npy')
        ind.append(len(ener_arr))
        if (np.ma.size(ener_arr, axis = 0)<lr-1):
            pass
        else:       
            arrsave1.append(ener_arr[0,0])
            arrsave1.append(ener_arr[lr-2,0])
    arrsave1 = np.reshape(arrsave1, (int((len(arrsave1))/2),2), order = 'C')
    return arrsave1, ind#, arrsave2

