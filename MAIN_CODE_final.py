

"""This code will solve the Berger and Seltzer stopping power ode for singular time 
but with an artbitrary pellet potential which linearly increases from the pellet to 
the cloud edge.

A potential field can be applied to account for additional electrostatic slowing of particles
and changes in this potential to plot how the fraction of the initial distribution changes 
with this pellet potential """ 
import numpy as np
from gen_var import rp, rc, t_start,lt, lr, dr, lp, pel_pot, cloud_pot, style , many_start, t_end, inc, p_inc, phi_sheath
import os
import stopblock #function to bethe stop electrons
import MB_calc # function to determine EEDF
import pot_sum # function to add potentials from electrons and pellet at matching points
from electron import e_dist, ener_res, e_bar, KE_top, KE_bot, RME, M_fac #importing essential variables and functions
import electron as part
import bethe_analysis_functions as baf # funcitons to analyse the data
import sys #for exiting the code and printing an "error message"
import discret_mat
import iterative_sol as SOR 
import elec_transport_2 as e_trans #THIS NEEDS TO BE PUT IN YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#look closer at this transport function
import common_interpolator as com_int 
#import cloud_mover 



particle = 'electron' # clarifying particle type

e_mid, e_bins, MB_norm = part.dist_calc(e_dist, ener_res, e_bar) #calculating mid-point energy of bins, bins, the normalised distribution
le = len(e_mid)
ener = np.arange(KE_bot, KE_top, ener_res) # establishing energy array
MB = MB_calc.MB(ener, e_bar) 

#############################################

i = t_start #index of time-array

#############################################

if (i > lt-1):
    print('Index for time array is too large. Select a time-index that \
        fits within length of t array')
    sys.exit()
else:
    pass

r = np.arange(0,rc[-1], dr) 

lr = len(r)

"Must now establish how the simulation will proceed with choice of style"
if style == 'neutral':
    savedir_raw = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/one_iteration/raw_outputs/'
    savedir_an = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/one_iteration/analysed_outputs/'
elif style =='self_charge':
    savedir_raw = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/one_iteration_phic/raw_outputs/'
    savedir_an = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/one_iteration_phic/analysed_outputs/'
elif style == 'self_sheath_charge':
    savedir_raw = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/many_iteration/raw_outputs/'
    savedir_an = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/many_iteration/analysed_ouputs/'
else:
    print ('Must select an appropriate simulation style - see gen_var for details')
    sys.exit()

"""Delete all files in the sub directory"""

filelist_raw = [ f for f in os.listdir(savedir_raw) if f.endswith(".npy") ]
for f in filelist_raw:
    os.remove(os.path.join(savedir_raw, f))

filelist_an = [ f for f in os.listdir(savedir_an) if f.endswith(".txt") ]
for f in filelist_an: 
    os.remove(os.path.join(savedir_an, f))

"Establish the new arrays"
pot = np.zeros(lr)
flux_arr = []
ener_flux_arr = []
lifetime_arr = []
acc_elec_dens = np.zeros(len(r)) # initial electron density - zero everywhere
pot /= RME*M_fac # normalise electric potential for stopblock calc

rp_crit = 0.001 # can get put in gen_var later

low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
if style == 'self_sheath_charge':
    pot = np.linspace(-phi_sheath, 0.0, r[low:up])
while rp[i] > rp_crit:
        low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
        up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
        #r_internal = np.arange(0, rc[i] - rp[i] , dr)
        r_internal = r[low:up]
        stopblock.stopblock_phi(e_mid, r_internal,i, particle, pot)

        term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r_internal, i, le, savedir_raw)
        stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r_internal,i,len(e_mid), savedir_raw)
        faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,p, r_internal[0])
        ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,savedir_an, term_en)

        np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
        flux_arr.append((ret_flux_frac, pel_pot[p]))
        ener_flux_arr.append((ener_flux, pel_pot[p]))
        lifetime_arr.append((lifetime, pel_pot[p]))

        life = life + time_diff*delta_t

        elec_interp, ind_low, ind_up = com_int.common_interp(r, real_density[:,1], real_density[:,0]) # interpolating charge onto all gridpoints
        acc_elec_dens[ind_low:ind_up] = elec_interp[:] + acc_elec_dens[ind_low:ind_up]
        #now need to move the electrons
        acc_elec_dens, low_2, up_2 = e_trans.elec_mover(i,r,acc_elec_dens,r[ind_low:ind_up],time_diff)
        r_domain = r[low_2:up_2] #  shifted spatial range, rp2 to tc2
        if style == 'neutral':
            pass
        else:
            A = discret_mat.discret(r_domain) #r here needs to be from just in front of the pellet to the cloud at the new time.
            pot = SOR.SOR(A, pot[low_2:up_2], acc_elec_dens[low_2:up_2], r_domain)

        "Saving data for a singular Bethe calculation"

        np.savetxt(os.path.join(savedir_an, 'terminal_energy_pot_test_t'+str(i)+'pot'+str(pot[0]) +'.txt'), term_en)   
        np.savetxt(os.path.join(savedir_an, 'stop_point_pot_test_t' + str(i) +'pot'+str(pot[0])+'.txt'), stop_point, fmt = ('%f'))
        np.savetxt(os.path.join(savedir_an, 'density_pot_test_t' +str(i) +'pot'+str(pot[0])+'.txt'), faux_density)
        np.savetxt(os.path.join(savedir_an, 'real_density_pot_test_t'+str(i) +'pot'+str(pot[0])+'.txt'), real_density)



"Saving compilation of Bethe calculations with varying potentials"
np.savetxt(os.path.join(savedir_an, 'retarded_flux_pot_test_t'+str(i)+'.txt'), flux_arr)    
np.savetxt(os.path.join(savedir_an, 'retarded_ener_flux_pot_test'+str(i)+'.txt'), ener_flux_arr)
np.savetxt(os.path.join(savedir_an, 'pellet_lifetime_retarding_potential_pot_test' + str(i)+'.txt'), lifetime_arr)
print('success')