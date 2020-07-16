# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 16:23:56 2020

@author: Kyle
"""

import scipy.interpolate as spint
import numpy as np
def smoother(pos,bins,ind_up, ind_low):
    indices = []
    fac = 1.05
    for i in range(ind_low + 1, ind_up):
        if bins[i] > fac*bins[i+1]:
            indices = np.append(indices,i)
    for i in range(ind_up, len(pos)-1):
        if bins[i] > fac*bins[i+1]:
            indices = np.append(indices, i)
    bins_new = np.delete(bins,indices)
    pos_new = np.delete(pos, indices)

    f = spint.interp1d(pos_new, bins_new, kind = 'cubic')
    bins_new_mid = spint.pchip_interpolate(pos_new,bins_new,pos[ind_low+1:ind_up])
    bins_empty = np.zeros(ind_low)
    bin_pellet = bins[ind_low]
    bins_cloud = spint.pchip_interpolate(pos_new,bins_new,pos[ind_up:])
    bins_new = np.append(bins_empty,bin_pellet)
    bins_new = np.append(bins_new, bins_new_mid)
    bins_new = np.append(bins_new,bins_cloud)
    return bins_new