# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 19:19:36 2020

@author: Kyle
"""
import numpy as np 
import matplotlib.pyplot as plt 
import romberg as ro 
from gen_var import delta_t, vrms_e, dens_plas, r0, tf 

import os 

direc = os.getcwd() 

load_dir = os.path.join(direc, 'many_iteration','neutral','1000eV','analysed_outputs') + os.sep
file = np.genfromtxt(load_dir + 'outputs1.txt', delimiter = ',', dtype = 'str')

r = file[5:, 0]
dens = file[5:,1]

r = np.asarray([float(a) for a in r])
dens = np.asarray([float(a) for a in dens])

r_lim = 0.91*r0

print('REAL UNITS')
r *= r0
dens *= dens_plas

N = np.pi*dens_plas*r_lim**2 * delta_t * tf * vrms_e 

A = 4.0*np.pi*ro.romberg_samp(dens * r*r,r )
dens_e = dens*N/A

print('Particle number is ' + str(N))

check_real = 4.0*np.pi*ro.romberg_samp(dens_e*r*r,r)
print('Check on particle number is ' + str(check_real))

print('NON-DIMENSIONAL UNITS')

r/= r0
dens /= dens_plas 

Nn = (r_lim/r0)**2
An = ro.romberg_samp(dens*r*r,r)
dens_en = dens*Nn/An

print('Particle number is ' + str(Nn))

check_nd = ro.romberg_samp(dens_en*r*r,r)
print('Check on particle number is ' + str(check_nd))

print('SCALED UNITS')

Nhat = np.pi*r0**2 * delta_t *tf * dens_plas * vrms_e 

Ahat = r0**3 * dens_plas  * 4.0 * np.pi 

dens_ehat = dens_plas *Nhat/Ahat 

check_sca  = 4.0*np.pi*r0**3 * dens_ehat
print('Particle number is ' + str(Nhat))

print('Check on particle number is ' + str(check_sca))

print('FINAL CHECK')
print('Check on Densities: n = ~n * ^n')
print('n = ' + str(dens_e[15]))
print('~n * ^n = ' + str(dens_en[15] * dens_ehat))
print('Check on Numbers: N = ~N * ^N')
print('N_real = ' + str(check_real))
print('~N * ^N = ' + str(check_nd * check_sca))

print('-------------CYLINDER METHOD-----------------')

print('REAL UNITS')
r *= r0
dens *= dens_plas

N = 0.25*np.pi*dens_plas*r_lim**2 * delta_t * tf * vrms_e 

A = np.pi*r_lim*r_lim*ro.romberg_samp(dens,r )
dens_e = dens*N/A

print('Particle number is ' + str(N))

check_real = np.pi*r_lim*r_lim*ro.romberg_samp(dens_e,r)
print('Check on particle number is ' + str(check_real))

print('NON-DIMENSIONAL UNITS')

r/= r0
dens /= dens_plas 

Nn = (r_lim/r0)**2
An = (r_lim/r0)**2 *ro.romberg_samp(dens,r)
dens_en = dens*Nn/An

print('Particle number is ' + str(Nn))

check_nd = (r_lim/r0)**2*ro.romberg_samp(dens_en,r)
print('Check on particle number is ' + str(check_nd))

print('SCALED UNITS')

Nhat = 0.25*np.pi*r0**2 * delta_t *tf * dens_plas * vrms_e 

Ahat = r0**3 * dens_plas * np.pi 

dens_ehat = dens_plas *Nhat/Ahat 

check_sca  = np.pi*r0**3 * dens_ehat
print('Particle number is ' + str(Nhat))

print('Check on particle number is ' + str(check_sca))

print('FINAL CHECK')
print('Check on Densities: n = ~n * ^n')
print('n = ' + str(np.amax(dens_e[:])))
print('~n * ^n = ' + str(np.amax(dens_en[:]) * dens_ehat))
print('Check on Numbers: N = ~N * ^N')
print('N_real = ' + str(check_real))
print('~N * ^N = ' + str(check_nd * check_sca))