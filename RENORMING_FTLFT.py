# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 09:21:13 2020

@author: Kyle
"""
import numpy as np 
import romberg as ro 
from gen_var import N0_elec, A0, dens_plas, r0, delta_t, tf, r_grid, rc_hr
from gen_var import many_start, ne_sca

def renorm(i,dens):
    N = rc_hr[i]**2
    integral = ro.romberg_samp(dens*r_grid*r_grid,r_grid)
    #print('Ratio of two N/A is ' + str(N/integral))
    dens_new = dens*N/integral
    print('Particles in is ' + str(N*N0_elec))
    check = ro.romberg_samp(dens_new*r_grid**2,r_grid) * N0_elec
    print('Contained charge = ' + str(check))
    check_again = 4.0*np.pi*ro.romberg_samp(dens_new * ne_sca * (r_grid*r0)**2, r_grid*r0)
    print('In real units, contained charge is ' + str(check_again))
    dens_new *= ne_sca
    return dens_new

import os 

direc = os.getcwd() 

load_dir = os.path.join(direc, 'many_iteration','neutral','1000eV','analysed_outputs') + os.sep
file = np.genfromtxt(load_dir + 'outputs1.txt', delimiter = ',', dtype = 'str')

r = file[5:, 0]
dens = file[5:,1]

r = np.asarray([float(a) for a in r])
dens = np.asarray([float(a) for a in dens])

dens_out = renorm(many_start, dens)
