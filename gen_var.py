#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 15:32:45 2019

@author: kyle
"""

#STYLE here is selected so as to determine whether or not a single time is evaluated with no prior history,
# a single time is evaluated with a prescribed potential, or whether time is permitted to advance with 
#electrostatic history either by a forced potential or the previously deposited charge distribution

#STYLE can be selected to be either 'once', 'once_charge', 'many'
style = 'many'

import numpy as np
import math as mt
import os
import stop_calc_rp_rc
#import RK4 
from electron import RME, M_fac, e_bar, e_dist, ener_res
import electron
import pellet_sheath as ps
#import rc_solve
import find_nearest as fn
import final_sheath_func as fsf
import scipy.interpolate as spint


direc = os.getcwd()
"""#TEMPORAL TERMS"""
tf = 1*10**(-3) #Pellet Lifetime - 1ms
tlow = 0.00
t_low_hr = np.copy(tlow) #rk45 solver will not work otherwise
tupper = 1.0
t_upper_hr = 0.2
dt = 0.002
dt_hr = 1e-6
t = np.arange(tlow, tupper, dt) # SETTING TIME GRID
t_hr = np.arange(t_low_hr,t_upper_hr,dt_hr)
lt = len(t)
lr = lt
t_start = 1
t_end = 500
inc = 1
many_start = 20
delta_t = 5*10**(-5) # normalised quantity - needed for number flux calculation
life = 0.0

""""#SPATIAL TERMS"""

r0 = 10**(-3) #Initial pellet radius in m
r0cgs = r0*100.0 #Initial pellet size in cm for normalisation in CSDA 
rp_crit = r0/(10**3) # normalised to initial pellet radius in m
rp = stop_calc_rp_rc.calc_rp(t) #CALCULATES PELLET RADIUS AT TIME POINTS IN NORMALISED UNITS 
#rp_hr = stop_calc_rp_rc.calc_rp(t_hr)
rc = stop_calc_rp_rc.calc_rc(t) #CALCULATES CLOUD RADIUS AT TIME POINTS IN NORMALISED UNITS 
#rc_hr = stop_calc_rp_rc.calc_rc(t_hr)
#rc_hr = rc_solve.rc_sol(t_hr)
rpd = np.zeros(len(t))
rpd_hr = np.zeros(len(t_hr))
for i in range(0, len(t)):
    rpd[i] = -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1)))))*(mt.tanh(1 - t[i])/(np.sqrt(np.log(mt.cosh(1 - t[i])))))

#for i in range(0, len(t_hr)):
#    rpd_hr[i] = -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1)))))*(mt.tanh(1 - t_hr[i])/(np.sqrt(np.log(mt.cosh(1 - t_hr[i])))))

rp_hr = spint.pchip_interpolate(t,rp,t_hr)
rc_hr = spint.pchip_interpolate(t,rc,t_hr)
rpd_hr = spint.pchip_interpolate(t,rpd,t_hr)


dr = 5.0*10**-1 #Spatial resolution 
#r_grid = np.linspace(0.0,rc[-1], rgl , endpoint = 'true')

n = 5
while dr > 5e-3:
    n += 1
    rgl = 2**n + 1
    r_grid = np.linspace(rp_hr[-1] , rc_hr[-1], num = rgl)
    dr = r_grid[1] - r_grid[0]
r_grid = np.arange(rp_hr[-1], rc_hr[-1], dr)
rgl = len(r_grid) # r grid lengths
"""#BERGER AND SELTZER/BETHE STOPPING POWER TERMS"""
zovera = 0.5 #Z/A
I = 19.8 #Mean excitation energy in eV, could be 19.8 from the molecular gaseous term from BERGER AND SELTZER 1984
c1= 0.153536 #Some constant which carries units of MeVcm^2/g
solid_dens = 0.086 #Solid density of Hydrogenin g/cm^3 
masscgs = solid_dens*r0cgs**3 #Units mass in g
delta = 0.0 # DEnsity correction, zero for H target

r_start = np.asarray([0.0])
eps = 0.01 #DEnsity contrast term - dimensionless
trunc_fac = 1.01*np.sqrt(2.0/np.exp(1))

"""NUMERICAL QUANTITIES"""


le = 500 #Number of energy points

m_e = 9.10938356*10**(-31) # in kg
m_p = 1.6726219*10**(-27) # in kg
"""ELECTROSTATIC TERMS"""

e = 1.6022*10**(-19) # elementary charge in Coulombs
epsilon0 = 8.8542*10**(-12) # permittivity of free space in F/m
phi_plas = 1.0*10**3 #plasma potential in V
phi_plas = e_bar
phi_p = -phi_plas*np.log(np.sqrt(2.0*m_p/(2.0*np.pi*m_e))) #Floating sheath potential in H2 plasma but need to look at a text book for this

p_diff = 200.0 # following terms needed for setting up a "pseudo-sheath" or potential across the cloud
pel_pot = -np.arange(0.0, 3000 + p_diff, p_diff)
cloud_pot = 0.0
lp = len(pel_pot)
p_inc = 1
sig = [0.5,0.6,0.7,0.8,0.9,1.0,1.25]
#sig = np.arange(0.40,1.50,0.05)

"SHEATH SOLUTION"
#s,sheath_pot, sheath_x = ps.sheath_solver()
s, sheath_pot, sheath_x = fsf.sheath()
s = np.abs(s)
s /=r0
many_start = fn.find_nearest(rc_hr - rp_hr,  s)

delta_t_new = s*r0*np.sqrt(m_e/(2.0*e*e_bar))
delta_t_new /= tf
#many_end = fn.find_nearest(rc_hr - rp_hr,3.0)

"""BACKGROUND PLASMA TERMS""" 
phi_plas = np.copy(e_bar) # in Volts
dens_plas = 10**19 # in m^-3
vrms_e = np.sqrt(2.0*e*phi_plas/m_e) #in ms^-1
F0 = 0.25*dens_plas*vrms_e # flux of particles per square metre per second


"""PELLET TERMS"""
bond_energy = 0.00527 #in eV
pel_dens_numb = solid_dens*10**(-3) # in kg/cm^3
pel_dens_numb *= 10**6
pel_dens_numb /= m_p #pellet number density in m^-3

"CROSS SECTIONS"
h2_ion_rate = np.loadtxt('h2_ion_rate.txt')
cs_ion_neut = 1e-20

"DIFFUSION TERMS" 
#s = 2.968 #Lennard Jones Average Collision Diameter in angstroms
omega = 0.1 # not exact but works for now, may requre interpolation to fix
A_con = 1.858*10**(-3) # in complicated units
m_h2 = 2.0*10**(-3) # molar mass for H2 approximately
m_h = 1.0*10**(-3) # molar mass of H approximately
r_rate = np.loadtxt(direc + os.sep + 'h2_ion_rate.txt')
mfp = (cs_ion_neut*pel_dens_numb)**-1
v_ion = np.sqrt(2*I*trunc_fac*e/(2.0*m_p))
mft = mfp/v_ion

diff_ion = mfp*mfp/mft
ratio = 10

"SCALED QUANTITIES FOR DIFFUSION"
diff_ion_sca = np.copy(diff_ion)
diff_term_sca = (r0/(np.sqrt(np.pi*diff_ion_sca*tf)))*np.exp(-r0**2/(diff_ion_sca*tf))
ion_flux_sca = diff_ion_sca*dens_plas/r0
ion_dens_renorm = diff_ion_sca*delta_t/(ratio*r0*r0)

"EEDF TERMS"
e_mid, e_bins, MB_norm = electron.dist_calc(e_dist, ener_res, e_bar)

"SCALED QUANTITIES FOR ABLATION" 

no_flux_sca = 0.25*dens_plas*vrms_e 
e_flux_sca = no_flux_sca*RME*M_fac 
ener_ablat_sca = e_flux_sca*r0**2 * tf *delta_t
N_ablat_sca = ener_ablat_sca/bond_energy
N_0_sca = (4.0/3.0)*np.pi*r0**3 * pel_dens_numb
mfp_ion = (cs_ion_neut*pel_dens_numb*eps)**-1
mfp_ener = np.abs(sheath_pot[0] - spint.pchip_interpolate(sheath_x,sheath_pot, mfp_ion) )

"scaling test"
Np_now = rp[20]**3 *(N_0_sca/N_ablat_sca)
Np = Np_now * N_ablat_sca