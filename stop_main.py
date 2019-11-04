#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 14:17:13 2019

@author: kyle
"""

#This new code will use modules to determine the stopping of a particle, electron or proton,
#in our ablation cloud using the modified Bethe stopping power form as determined by 
#Seltzer and Berger. This is no faster than other versions but it will be self-consistent
#and will only require this singular file to run and will call upon other modules where
#necessary. 

import os
import numpy as np
from gen_var import * # imports general variables
import stop_analysis_particle_density # calculates stops fraction of MB
import stop_analysis_stop_point # calculates position of stopped electrons
import stop_analysis_term_ener # calculates impacting energy of particles from distribution
import phi_calc # calculates the potential 
import MB_calc #calculates the maxwell boltzmann distribution
#Check for particle type and then choose what formula(modules) to apply
particle = 'electron'

if (particle =='electron'):
    from electron import *
    import electron as part

elif (particle == 'proton'):
    from proton import *
    import proton as part
else:
    print ('Particle type has no data. Please create particle data and re-run code') 
    

#import stop_func
import stopblock # code to solve dE/dx for initial energy

e_mid, e_bins, MB_norm = part.dist_calc(e_dist, ener_res, e_bar) #output bin densities, energy midpoints and normalised initial distribution

print(len(MB_norm))

ener = np.arange(100.0, 20000.0, ener_res)
MB = MB_calc.MB(ener, 1000.0)

print('Length of e_mid is' + str(len(e_mid)))
print('Length of initial distribution is' + str(le))
for i in range(1, lt):  
    print(i)
    r = np.arange(0, rc[i*int(len(rc)/lt)] - rp[i*int(len(rc)/lt)],dr)  #array of zero to pellet surface through cloud
    r = np.arange(0, rc[i] - rp[i], dr) #CHANGE MADE 30/1/2019
    """if (i==1):
       phi = np.zeros(len(r))
    else:
       pass"""
    stopblock.stopblock(e_mid,r,i, particle)

    term_en, ind = stop_analysis_term_ener.term_energy(particle, r, i, le)
    stop_point = stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid))
    faux_density = stop_analysis_particle_den sity.particle_density(stop_point,i, le, e_bins, particle)
    

    np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_terminal_ener_1keV_fs_t'+str(t) +'.txt'), term_en)  
    np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_deposition_cloud_1keV_fs_t' + str(i) +'.txt'), stop_point, fmt = ('%f')) 
    np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_density_cloud_1keV_fs_t' +str(t) +'.txt'), faux_density)
    #phi = phi_calc.phi_calc(i, particle, faux_density)
    
        
