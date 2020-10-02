# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:01:32 2020

@author: Kyle
"""
import numpy as np
def sorter(r,dens):
    r_sorted = np.sort(r)
    dens_sorted = np.zeros(len(r))
    for i in range(0, len(r)):
        index = np.argmin(np.abs(r_sorted[i] - r[:]))
        dens_sorted[i] = dens[index] 
    return r_sorted, dens_sorted
    