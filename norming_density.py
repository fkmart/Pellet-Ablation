# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 09:15:11 2020

@author: Kyle
"""

def norm(faux_density, r):
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

def norm_sort(faux_density, r):
    import numpy as np
    import scipy.interpolate as spint
    import sorter
    pos = faux_density[:-1,-1]
    bins = faux_density[:-1,0]
    pellet_bin = faux_density[-1,0]
    total = np.sum(bins)
    r_sort, dens_sort = sorter.sorter(pos,bins)
    
    #Remove any points that are co-spatial
    
    r_u = np.unique(r_sort) 
    dens_u = np.zeros(len(r_u))
    
    for a in range(0, len(dens_u)):
        for s in range(0, len(r_sort)):
            if r_u[a] == r_sort[s]: 
                dens_u[a] += dens_sort[s]
    dens_interp = spint.pchip_interpolate(r_u,dens_u,r[1:])
    dens_interp = dens_interp*(total/(np.sum(dens_interp)))
    dens_interp = np.append(pellet_bin,dens_interp)
    return dens_interp
    