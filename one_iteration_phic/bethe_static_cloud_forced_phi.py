import numpy as np
from gen_var import t_static, rp, rc, dr, dt, t, cloud_pot
from gen_var import pel_pot as pot
import os

"""This code will solve the Berger and Seltzer stopping power ode for singular time 
but with an artbitrary pellet potential which linearly increases from the pellet to 
the cloud edge.

A potential field can be applied to account for additional electrostatic slowing of particles
and changes in this potential to plot how the fraction of the initial distribution changes 
with this pellet potential """ 

import stopblock #function to bethe stop electrons
import MB_calc # function to determine EEDF

particle = 'electron' # clarifying particle type

from electron import e_dist, ener_res, e_bar, KE_top, KE_bot, RME, M_fac #importing essential variables and functions
import electron as part
import bethe_analysis_functions as baf

"""Delete all files in the sub directory"""

savedir_raw = os.path.abspath('/home/kyle/Documents/Python_Code/stop_code/one_iteration_phic/raw_outputs/')
savedir_analysed = '/home/kyle/Documents/Python_Code/stop_code/one_iteration_phic/analysed_outputs'
filelist = [ f for f in os.listdir(savedir_raw) if f.endswith(".npy") ]
for f in filelist:
    os.remove(os.path.join(savedir_raw, f))

filelist_analysed =  [ f for f in os.listdir(savedir_analysed) if f.endswith(".txt") ]
for f in filelist_analysed:
    os.remove(os.path.join(savedir_analysed, f))

e_mid, e_bins, MB_norm = part.dist_calc(e_dist, ener_res, e_bar) #calculating mid-point energy of bins, bins etc

ener = np.arange(KE_bot, KE_top, ener_res) # establishing energy array
MB = MB_calc.MB(ener, e_bar) 

flux_arr = []
ener_flux_arr = []
lifetime_arr = []

i = t_static

r = np.arange(0,rc[i] - rp[i], dr)



for p in range(0, len(pot)):
    print(pot[p])
    lr = len(r)
    phi = np.linspace(cloud_pot, pot[p], lr)
    phi[:] /=RME*M_fac
    stopblock.stopblock_phi(e_mid, r,i, particle,phi, savedir_raw)

#Analysis codes 

    term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, len(e_mid), savedir_raw)
    stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), savedir_raw)

    np.savetxt(savedir_analysed +'/terminal_energy_t'+str(i)+'pot' + str(pot[p]) +'.txt', term_en)   
    np.savetxt(savedir_analysed +  '/stop_point_t' + str(i)+ 'pot' + str(pot[p]) + '.txt', stop_point, fmt = ('%f'))

    faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle)
    ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(t_static, term_en)
    

    np.savetxt(savedir_analysed +  '/density_t' +str(i)+ 'pot' + str(pot[p]) + '.txt', faux_density)
    np.savetxt(savedir_analysed + '/real_density_t'+str(i)+ 'pot' + str(pot[p]) + '.txt', real_density)
    np.append(flux_arr, (ret_flux_frac, i*dt/(len(t)))) #Append potential dependant arrays with new quantities
    flux_arr.append((ret_flux_frac, pot[p]))
    ener_flux_arr.append((ener_flux, pot[p]))
    lifetime_arr.append((lifetime, pot[p]))
#stop time loop here
"Saving compilation of Bethe calculations with varying potentials"
np.savetxt(savedir_analysed +  '/retarded_flux_t'+str(i) + '.txt', flux_arr)
np.savetxt(savedir_analysed + '/retarded_ener_flux'+str(i)+ '.txt', ener_flux_arr)
np.savetxt(savedir_analysed +  '/pellet_lifetime_retarding_potential' + str(i) + '.txt', lifetime_arr)
print('success')