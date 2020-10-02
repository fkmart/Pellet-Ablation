# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:32:04 2020

@author: Kyle
"""
import numpy as np 
from gen_var import cs_ion_neut, v_ion, dr, m_p, e, RME, M_fac, r0

def crossing(mfp, ind_low_r, ind_up_r,phi,delta_t_ion):
    time = 0.0  # add to this until time equals delta_t_ion
    k = 0
    t_cross = 0.0
    while phi[ind_low_r + k] < 0.0:
        k +=1
    j = ind_low_r + 1
    k += ind_low_r
    while j < k and time < delta_t_ion:
        acc = np.abs((phi[j-1] - phi[j])/(dr*2.0*m_p)) * e * RME*M_fac
        t_cross = -(v_ion - np.sqrt(v_ion**2 + 2.0*acc*mfp[j]))/(acc)
        t_cross *= r0*dr/mfp[j]
        time += t_cross
        j +=1
    return j 
    