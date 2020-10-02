# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:44:45 2020

@author: Kyle
"""

import numpy as np 
import scipy.interpolate as spint
def cs_calc(ener):
    file = np.loadtxt('phelps_ion_neut.txt')
    energies = file[:,0]
    cs = file[:,1]
    cs_ion_neut = spint.pchip_interpolate(energies,cs, ener)
    return cs_ion_neut 
