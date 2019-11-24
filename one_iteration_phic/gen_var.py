#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 15:32:45 2019

@author: kyle
"""

import numpy as np
#from gen_var import *
#rp = np.load('rp_500.npy') 
#rc = np.load('rc500.npy') 

#from stop_calc_rp_rc import *
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
t_static = 80

""""#SPATIAL TERMS"""

r0 = 10**(-3) #Initial pellet radius in mm
r0cgs = 10**-1 #Initial pellet size in cm for normalisation
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

e = 1.6*10**(-19) # elementary charge in Coulombs
epsilon0 = 8.854*10**(-12) # permittivity of free space in F/m
EF0 = -e/(4*np.pi*epsilon0*r0**2) # electric field at critical points for singular charge
phi_plas = 1.0*10**3 #plasma potential in eV
phi_p = -2.83*phi_plas #Floating sheath potential in H2 plasma but need to look at a text book for this

pel_pot = -np.arange(0.0, 3020.0, 20.0)
cloud_pot = 0.0
lp = len(pel_pot)
sig = [0.70,0.75,0.80,0.85,0.90,0.95,1.00]

"""Background plasma terms""" 
phi_plas = 10**3 
m_e = 9.11*10**(-31)
m_p = 1.67*10**(-27)
dens_plas = 10**20 
vrms_e = np.sqrt(2.0*e*phi_plas/m_e)
F0 = 0.25*dens_plas*vrms_e # flux of particles per square metre per second


"""Pellet terms"""
bond_energy = 0.00527 #in eV
pel_dens_numb = solid_dens*10**(-3) # in kg/cm^3
pel_dens_numb *= 10**6
pel_dens_numb /= m_p 
