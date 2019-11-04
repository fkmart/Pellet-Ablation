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

pot_test = 'false'

particle = 'electron' # clarifying particle type

from electron import e_dist, ener_res, e_bar, KE_top, KE_bot, RME, M_fac #importing essential variables and functions
import electron as part
import bethe_analysis_functions as baf

"""Delete all files in the sub directory"""

mydir = './static_outputs_phi'
filelist = [ f for f in os.listdir(mydir) if f.endswith(".npy") ]
for f in filelist:
    os.remove(os.path.join(mydir, f))

e_mid, e_bins, MB_norm = part.dist_calc(e_dist, ener_res, e_bar) #calculating mid-point energy of bins, bins etc

ener = np.arange(KE_bot, KE_top, ener_res) # establishing energy array
MB = MB_calc.MB(ener, e_bar) 

#############################################

i = t_static #index of time-array

#############################################

if (i > lt-1):
    print('Index for time array is too large. Select a time-index that \
        fits within length of t array')
else:
    pass

r = np.arange(0,rc[i] - rp[i], dr)

lr = len(r)
"Now set up the potential field across the cloud" 

flux_arr = []
ener_flux_arr = []
lifetime_arr = []

if (pot_test == 'true') :
    elec_pot = np.loadtxt(mydir+'/elec_pot_t_static.txt')
    r_elec = elec_pot[:,1]
    elec_pot = elec_pot[:,0]
else:
    pass

for p in range(0, lp):
    pot = np.linspace(pel_pot[p], cloud_pot, lr)
    if (pot_test =='true'):
        pot = pot_sum.pot_sum(elec_pot, pot, r_elec, r)
    else:
        pass
    pot /= RME*M_fac
    stopblock.stopblock_phi(e_mid, r,i, particle, pot)

    #Analysis codes 

    term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, mydir)
    stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), mydir)
   
    
    
    "Saving data for a singular Bethe calculation"
    if (pot_test =='true'):
        np.savetxt(os.path.join(mydir, 'terminal_energy_pot_test_t'+str(i)+'pot'+str(pel_pot[p]) +'.txt'), term_en)   
        np.savetxt(os.path.join(mydir, 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'.txt'), stop_point, fmt = ('%f'))
    else:
        np.savetxt(os.path.join(mydir, 'terminal_energy_t'+str(i)+'pot'+str(pel_pot[p]) +'.txt'), term_en)   
        np.savetxt(os.path.join(mydir, 'stop_point_t' + str(i) +'pot'+str(pel_pot[p])+'.txt'), stop_point, fmt = ('%f'))
    faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,p)
    ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(t_static,mydir, term_en)
    
    if (pot_test =='true'):
        np.savetxt(os.path.join(mydir, 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'.txt'), faux_density)
        np.savetxt(os.path.join(mydir, 'real_density_pot_test_t'+str(i) +'pot'+str(pel_pot[p])+'.txt'), real_density)
    else:
        
        np.savetxt(os.path.join(mydir, 'density_t' +str(i) +'pot'+str(pel_pot[p])+'.txt'), faux_density)
        np.savetxt(os.path.join(mydir, 'real_density_t'+str(i) +'pot'+str(pel_pot[p])+'.txt'), real_density)
    np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
    flux_arr.append((ret_flux_frac, pel_pot[p]))
    ener_flux_arr.append((ener_flux, pel_pot[p]))
    lifetime_arr.append((lifetime, pel_pot[p]))
"Saving compilation of Bethe calculations with varying potentials"
if (pot_test =='true'):
    np.savetxt(os.path.join(mydir, 'retarded_flux_pot_test_t'+str(i)+'.txt'), flux_arr)
    np.savetxt(os.path.join(mydir, 'retarded_ener_flux_pot_test'+str(i)+'.txt'), ener_flux_arr)
    np.savetxt(os.path.join(mydir, 'pellet_lifetime_retarding_potential_pot_test' + str(i)+'.txt'), lifetime_arr)
else:
    np.savetxt(os.path.join(mydir, 'retarded_flux_t'+str(i)+'.txt'), flux_arr)
    np.savetxt(os.path.join(mydir, 'retarded_ener_flux'+str(i)+'.txt'), ener_flux_arr)
    np.savetxt(os.path.join(mydir, 'pellet_lifetime_retarding_potential' + str(i)+'.txt'), lifetime_arr)
print('success')