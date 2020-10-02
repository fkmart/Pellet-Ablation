# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 22:46:56 2020

@author: Kyle
"""
import numpy as np
import romberg as ro 
from gen_var import r_grid, rc_hr, N0_elec, dens_plas
def renorm(i,dens,N_tot):
    integral = 4.0*np.pi*ro.romberg_samp(dens*r_grid*r_grid,r_grid)
    dens_new = dens*frac*N0_elec/integral
    return dens_new