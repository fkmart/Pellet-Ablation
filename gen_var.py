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
style = 'once'

import numpy as np


import stop_calc_rp_rc
import RK4 


"""#TEMPORAL TERMS"""
tf = 1*10**(-3) #Pellet LIfetime - 1ms
tlow = 0.0

tupper = 1.0

dt = 0.002
t = np.arange(tlow, tupper, dt) # SETTING TIME GRID
lt = len(t)
lr = lt
t_start = 50
t_end = 500
inc = 100
many_start = 3

""""#SPATIAL TERMS"""

r0 = 10**(-3) #Initial pellet radius in mm
r0cgs = 10**-1 #Initial pellet size in cm for normalisation in CSDA 
rp_crit = 0.001 # normalised to initial pellet radius in mm 
rp = stop_calc_rp_rc.calc_rp(t) #CALCULATES PELLET RADIUS AT TIME POINTS IN NORMALISED UNITS 
rc = stop_calc_rp_rc.calc_rc(t) #CALCULATES CLOUD RADIUS AT TIME POINTS IN NORMALISED UNITS 

"""#BERGER AND SELTZER/NETHE STOPPING POWER TERMS"""
zovera = 0.5 #Z/A
I = 22.3 #Mean excitation energy in eV
c1= 0.153536 #Some constant which carries units of MeVcm^2/g
solid_dens = 0.086 #Solid density of Hydrogenin g/cm^3 
masscgs = solid_dens*r0cgs**3 #Units mass in g
delta = 0.0 # DEnsity correction, zero for H target

r_start = np.asarray([0.0])
eps = 0.01 #DEnsity contrast term - dimensionless

"""NUMERICAL QUANTITIES"""

dr = 1.0*10**-2 #Spatial resolution 
le = 500 #Number of energy points

"""#electrostatic terms"""

e = 1.6022*10**(-19) # elementary charge in Coulombs
epsilon0 = 8.8542*10**(-12) # permittivity of free space in F/m
EF0 = -e/(4*np.pi*epsilon0*r0**2) # electric field at critical points for singular charge
phi_plas = 1.0*10**3 #plasma potential in eV
phi_p = -2.83*phi_plas #Floating sheath potential in H2 plasma but need to look at a text book for this

p_diff = 1000.0 # following terms needed for setting up a "pseudo-sheath" or potential across the cloud
pel_pot = -np.arange(0.0, 3000 + p_diff, p_diff)
cloud_pot = 0.0
lp = len(pel_pot)
p_inc = 1

"""Background plasma terms""" 
phi_plas = 10**3 # in Volts
m_e = 9.10938356*10**(-31) # in kg
m_p = 1.6726219*10**(-27) # in kg
dens_plas = 10**20 # in m^-3
vrms_e = np.sqrt(2.0*e*phi_plas/m_e) #in ms^-1
F0 = 0.25*dens_plas*vrms_e # flux of particles per square metre per second


"""Pellet terms"""
bond_energy = 0.00527 #in eV
pel_dens_numb = solid_dens*10**(-3) # in kg/cm^3
pel_dens_numb *= 10**6
pel_dens_numb /= m_p #pellet number density in m^-3
