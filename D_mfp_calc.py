# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 12:12:51 2020

@author: Kyle
"""
import numpy as np 
from gen_var import cs_ion_neut, v_ion, pel_dens_numb,eps, r_grid, rp_hr

def mfp_calc(ind_low, ind_up, i):
    dens_neut = pel_dens_numb*eps*(1.0 + rp_hr[i]**2)/(1.0 + r_grid[ind_low:ind_up]**2)
    mfp_ion = (dens_neut*cs_ion_neut)**(-1)
    return mfp_ion

def diff_calc(mfp):
    diff_ion = mfp* v_ion
    return diff_ion 
