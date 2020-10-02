import numpy as np
from gen_var import * 
import os

"""This code will solve the Berger and Seltzer stopping power ode for singular time 
but with an artbitrary pellet potential which linearly increases from the pellet to 
the cloud edge.

A potential field can be applied to account for additional electrostatic slowing of particles
and changes in this potential to plot how the fraction of the initial distribution changes 
with this pellet potential """ 

import stopblock #function to bethe stop electrons
import MB_calc # function to determine EEDF
import pot_sum # function to add potentials from electrons and pellet at matching points

particle = 'electron' # clarifying particle type

from electron import e_dist, ener_res, e_bar, KE_top, KE_bot, RME, M_fac #importing essential variables and functions
import electron as part
import bethe_analysis_functions as baf

"""Delete all files in the sub directory"""

mydir = './multi_static_outputs'
myotherdir = './raw_multi_static_outputs'
filelist = [ f for f in os.listdir(myotherdir) if f.endswith(".npy") ]
for f in filelist:
    os.remove(os.path.join(myotherdir, f))

e_mid, e_bins, MB_norm = part.dist_calc(e_dist, ener_res, e_bar) #calculating mid-point energy of bins, bins etc

ener = np.arange(KE_bot, KE_top, ener_res) # establishing energy array
MB = MB_calc.MB(ener, e_bar) 

#############################################

i = t_static #index of time-array

#############################################

"Now set up the potential field across the cloud" 

flux_arr = []
ener_flux_arr = []
lifetime_arr = []

inc = 5 
#######################################
#NEED AN ARRAY HERE TO WORK THROUGH DIFFERENT TIMES
#######################################
for i in range(inc, lt, inc):
    r = np.arange(0,rc[i] - rp[i], dr)
    lr = len(r)
    stopblock.stopblock(e_mid, r,i, particle)

    #Analysis codes 

    term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, myotherdir)
    stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), myotherdir)
   
    "Saving data for a singular Bethe calculation"
    np.savetxt(os.path.join(mydir, 'terminal_energy_t{:03d}.txt'.format(i)   ), term_en)
    np.savetxt(os.path.join(mydir, 'stop_point_t{:03d}.txt'.format(i)), stop_point)
    faux_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,myotherdir)
    np.savetxt(os.path.join(mydir, 'density_t{:03d}.txt'.format(i)), faux_density)
    #Need to re-arrange faux_density array for use in normalising function 

    a1 = np.flip(faux_density[:,0], axis = 0)
    a2 = np.flip(faux_density[:,1], axis = 0)
    real_density = baf.stop_analysis_norming_elec_density.renorm_dens(a1, a2, e_bins, mydir, stop_point)

    #ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(t_static,myotherdir, term_en)
    np.savetxt(os.path.join(mydir, 'real_density_t{:03d}.txt'.format(i)), real_density)
    """
    #MORE WORK NEEDED TO FIX THIS
    np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
    flux_arr.append((ret_flux_frac, pel_pot[p]))
    ener_flux_arr.append((ener_flux, pel_pot[p]))
    lifetime_arr.append((lifetime, pel_pot[p]))
"Saving compilation of Bethe calculations with varying potentials"
    np.savetxt(os.path.join(mydir, 'retarded_flux_t'+str(i)+'.txt'), flux_arr)
    np.savetxt(os.path.join(mydir, 'retarded_ener_flux'+str(i)+'.txt'), ener_flux_arr)
    np.savetxt(os.path.join(mydir, 'pellet_lifetime_retarding_potential' + str(i)+'.txt'), lifetime_arr)"""
print('success')