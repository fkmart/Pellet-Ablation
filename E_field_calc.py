#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 14:49:12 2019

@author: kyle
"""

import numpy as np
import scipy.interpolate as spiter
import scipy.integrate as spiteg
import sphere_integrand
from gen_var import *

#STILL TO BE ADDED FOR CONSISENCY IS THE TOTAL CHARGE ON THE PELLET DUE TO SHEATH

E_perm = 8.854*10**(-12)

def E_field(cloud_dens, r, rp, i):
    def find_nearest(arr,val): #find nearest value in array
        array = np.asarray(arr)
        idx = np.abs(array - val).argmin()
        return idx
    elec_field = []
    elec_dens = cloud_dens[:,1] # taking deposited electron density from array
    f = spiter.interp1d(cloud_dens[:,0], elec_dens,  kind = 'cubic') #set up interpolation function
    ind_up = find_nearest(r,cloud_dens[0,0]) #index in r of first stopping position (cloud-plasma interface)
    ind_temp = find_nearest(cloud_dens[:,0], dr) #index in density profile of smallest step from pellet-cloud surface 
    ind_low = find_nearest(r, cloud_dens[ind_temp,0]) #index in r of value determined from above
    print(ind_low)
    print(ind_up)
    print(len(r))
    print('-----------------------------------')
    elec_dens = f(r[ind_low +2:ind_up]) # interpolate eletron density to r from deposition points
    #NEED TO DO THIS IN ONE DIMENSION FIRST. JUST CONSIDER A LINE TO START WITH
    #print('Total electron density is ' + str(spiteg.simps(elec_dens, r[ind_low+2:ind_up])))
    
    
    
    r = (r+rp[i])*10**-3
    
    for j in range(1, len(elec_dens)-1):
        #integ = 4*np.pi*(r[j]**2)*10**20 * np.flip(elec_dens[0:j], axis = 0) #this is the "3d" version. Symmetry should cancel theta and phi effects but may be needed for scaling
        integ_inside = np.flip(elec_dens[0:j], axis = 0)
        integ_outside = np.flip(elec_dens[j: len(elec_dens)-1], axis = 0)
        
        
        #print(len(integ))
        #x = spiteg.simps(integ_inside/E_perm, r[0:j]) - spiteg.simps(integ_outside/E_perm, r[j:len(elec_dens)-1])
        y = sphere_integrand.sph_integ_int(elec_dens[0:j], r[0:j])
        #print(x)
        elec_field = np.append(-e*elec_field/epsilon0,x) # convert EF into semi-real units
        elec_field = elec_field/EF0 # fully normalise
    return elec_field, elec_dens 

#TO REFLECT ACTUAL GEOMETRY
#CHECK REALITY OF NUMBERS AFTER THAT
