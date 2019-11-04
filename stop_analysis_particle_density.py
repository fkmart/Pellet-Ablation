#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:31:32 2019

@author: kyle
"""

import numpy as np 
import os
from gen_var import *
import stop_analysis_norming_elec_density 

def particle_density(arrsave2, t, le, ebins,particle,mydir,r_pel,r):
    e_dens = 0.0
    e_depo = []
    #posit = []
    arrsave3 = []
    arrsave4 = []
    for j in range(0,le -2):
         if arrsave2[j,1] == arrsave2[j+1,1]:
             e_dens = ebins[j] + e_dens

         elif arrsave2[j,1] != arrsave2[j+1,1]:
             e_dens = ebins[j] + e_dens
             e_depo.append(e_dens)
             arrsave3.append(np.asarray((arrsave2[j,1],e_dens)))
             e_dens = 0.0
    if arrsave2[-2,1] == arrsave2[-1,1]:
        pass
    else:
        arrsave3.append(np.asarray(arrsave2[-1,1]), ebins[-1])
    """if (e_dens ==0.0):
        arrsave3 = np.asarray(arrsave3)
    else:
        arrsave3 = np.asarray(e_dens)
        arrsave3 = np.append(arrsave3, arrsave2[j,1])
        arrsave3 = np.flip(arrsave3, axis = 0)"""
    arrsave3 = np.asarray(arrsave3)
    a1 = np.flip(arrsave3[:,0], axis = 0)
    a2 = np.flip(arrsave3[:,1],axis = 0)
    arrsave4 = stop_analysis_norming_elec_density.renorm_dens(a1[1:], a2[1:],ebins, arrsave2,t,r)
    #print(len(arrsave3))
   
    #np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_density_cloud_1keV_fs_t' +str(t) +'.txt'), arrsave3)
    return arrsave3, arrsave4
    