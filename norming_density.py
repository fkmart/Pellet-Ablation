# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 09:15:11 2020

@author: Kyle
"""

def norm(i,faux_density, r):
    #from gen_var import rp, rc
    import numpy as np 
    import scipy.interpolate as spint
    pos = faux_density[:-1,1]
    bins = faux_density[:-1,0]
    pellet_bin = faux_density[-1,0]
    total = np.sum(bins)
    g = spint.interp1d(pos,bins, kind = 'cubic', fill_value = 'extrapolate')
    #r = np.linspace(rp[i], rc[i], num = 513, endpoint = 'true')
    bins_interp = g(np.flip(r[:-1],axis = 0))
    pos_flipped = np.flip(pos,axis = 0)
    bins_flipped = np.flip(bins, axis = 0)
    bins_interp = spint.pchip_interpolate(pos_flipped,bins_flipped, r[1:])
    bins_interp = bins_interp*(total/np.sum(bins_interp))
    bins_final = np.append(pellet_bin,bins_interp)
    bins_final = np.flip(bins_final,axis = 0)
    return bins_final