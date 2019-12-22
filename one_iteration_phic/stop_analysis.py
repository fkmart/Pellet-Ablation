#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 12:35:28 2019

@author: kyle
"""

#New code to perform the analysis of the stopping code output.

import numpy as np

from gen_var import *
import os
#rp = np.load('rp_500.npy')
#rc = np.load('rc500.npy')

particle = 'electron'
print('Particle type is ' +str(particle))

if (particle == 'electron'):
    from electron import *
elif (particle == 'proton'):
    from proton import *
else:
    print('Particle is not known. Please create data file for particle and re-run')
    
le = len(np.arange(KE_bot, KE_top, ener_res))
import stop_analysis_term_ener
import stop_analysis_stop_point
import stop_analysis_particle_density
import E_field_calc

e_mid, e_bins, MB_norm = dist_calc(e_dist,ener_res, e_bar)

le_mid = len(e_mid)
cwd = os.getcwd()
#mydir = './Analysed_Outputs'
mydir = os.path.join(cwd, 'Analysed_Outputs') + os.sep

rp = stop_calc_rp_rc.calc_rp(t)
rc = stop_calc_rp_rc.calc_rc(t)
min_imp_ener = []
for t in range(1,lt):
    #print(t)
    r = np.arange(0, rc[t] - rp[t],dr)
    #print(r)
    #array,array2   = stop_analysis_term_ener.term_energy(particle,r,t,le_mid)
    
    #array3 = stop_analysis_particle_density.particle_density(array2,t,le_mid,e_bins, particle)

    term_en, ind = stop_analysis_term_ener.term_energy(particle, r, t, le, mydir)
    stop_point = stop_analysis_stop_point.stop_point(term_en,ind, particle, r,t,le_mid, mydir)
    faux_density = stop_analysis_particle_density.particle_density(stop_point,t, le, e_bins, particle)

    #e_field,dens_out = E_field_calc.E_field(array3, r, rp, t)

    #Still need to get the lowest impacting energies with time to show that shielding improves with time
    #It is a simple conclusion to draw, but effective in illustrating the decrease in particle numbers arriving at the 
    #pellet surface. Does not indicate the total amount of energy impacting on pellet. 

    if (particle=='electron'):
        print(array[0,0])
        min_imp_ener.append(array[0,0])
        min_imp_ener.append(t*dt)
    else:
        pass
if (particle =='electron'):
    array4 = np.reshape(min_imp_ener, (int(0.5*len(min_imp_ener)),2), 'C')
    np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_min_imp_ener_1keV_fs_t.txt'), array4)
else:
    pass