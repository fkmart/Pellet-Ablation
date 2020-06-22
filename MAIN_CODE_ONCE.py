"""This code will use an RK4 solver to solve the CSDA equation from Sletzer and Berger to determine the distribution
of electrons across the cloud due to inelastic collisions with neutral diatomic hydrogen. The distribution of charge
is then used to determine the self-field according to Poisson's equation by using an SOR method. The incident 
energy on the pellet surface then determines the """
import numpy as np
from gen_var import rp, rc, t_start,lt, lr, dr, lp, pel_pot, cloud_pot, style,r0, dens_plas,e,epsilon0
from gen_var import many_start, t_end, inc, p_inc, rp_crit, delta_t, lp , sig, r, n_r, r_grid, rgl
import os 
import stopblock #function to bethe stop electrons
import MB_calc # function to determine EEDF
from electron import e_dist, ener_res, e_bar, KE_top, KE_bot, RME, M_fac #importing essential variables and functions
import electron as part
import bethe_analysis_functions as baf # funcitons to analyse the data
import sys #for exiting the code and printing an "error message"
import discret_mat # makes discretisation matrix for SOR
import iterative_sol as SOR #SOR solver
import elec_transport_2 as e_trans #transport module to move electrons based on relative change in H2 densities
import common_interpolator as com_int #interpolater to move from one grid to another
import matplotlib.pyplot as plt # plotting module - not to be used in these simulations but present for testing
import gauss_test_pot # gaussian test potential to see effect of bump in potential
import gauss_centre_finder as gcf
from sig_calc import sig_calc
import grid_pusher as gp
import elec_transport_push 

particle = 'electron' # clarifying particle type
life = 0.0
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

neut_dens = np.zeros(n_r) # neutral density array of equal length, only some values will be non-zero
push_e_dens = np.zeros(rgl)
lr = len(r)

#Important to print the style being used here so you don't totally miss it
print('The scenario tested here is: ' + str(style))

"Must now establish how the simulation will proceed with choice of style"

pot = np.zeros(lr)
#pot = np.linspace(-2.8*10**-3, 0.0, num = lr, endpoint = 'true')
savedir_raw = os.path.join(os.getcwd(), 'one_iteration', 'raw_outputs') + os.sep
savedir_an = os.path.join(os.getcwd(), 'one_iteration', 'analysed_outputs') + os.sep

"""Delete all files in the sub directory"""
"""
filelist_raw = [ f for f in os.listdir(savedir_raw ) if f.endswith(".npy") ] # deletes the .npy files of raw data
for f in filelist_raw:
    os.remove(os.path.join(savedir_raw, f))

filelist_an = [ f for f in os.listdir(savedir_an) if f.endswith(".txt") ] # deletes the .txt files of analysed data
for f in filelist_an:
    os.remove(os.path.join(savedir_an, f))"""

"Establish the new arrays"
flux_arr = []
ener_flux_arr = []
lifetime_arr = []
acc_elec_dens = np.zeros(rgl) # initial electron density - zero everywhere
pot /= RME*M_fac # normalise electric potential for stopblock calc

p = 0
for i in range(10, 30, inc):
    print (i)
    save_raw_t = os.path.join(savedir_raw, 't_' + str(i)) + os.sep
    save_an_t = os.path.join(savedir_an, 't_' + str(i)) + os.sep
    if not os.path.exists(save_raw_t):
        os.mkdir(save_raw_t)
    if not os.path.exists(save_an_t):
        os.mkdir(save_an_t)
    low = next(p[0] for p in enumerate(r) if p[1] > rp[i]) # finds first index with r > r_p
    up = next(p[0] for p in enumerate(r) if p[1] > rc[i]) # finds first index with r > r_c
    r_internal = r[low:up] 
    neut_dens[low:up] = 0.01*((1.0 + rp[i]**2)/(1.0 + r[low:up]**2))
    r_grid = np.linspace(rp[i], rc[i], num = 512, endpoint = 'true')
    pot = np.zeros(len(r_grid))
    #pot[low:up] = np.linspace(-2.8e3,0.0,num = up-low, endpoint = 'true')
    #pot = np.linspace(-2.8e3,0.0, num = len(r_grid), endpoint = 'true')
    pot = pot/(RME*M_fac)
    stopblock.stopblock_phi_mod(e_mid, r,i, neut_dens, pot,save_raw_t) # change this line for new mods
    #stopblock.stopblock_phi_mod_rkf(e_mid, r_grid, i, pot, save_raw_t)
    term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, save_raw_t)
    stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), save_raw_t)
    #Now determine what index the particles reach the pellet
    ind = np.where(stop_point[:,1] < rp[i])
    ind = ind[0]
    faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,p,r)
    push_e_dens = gp.pusher(faux_density, r_grid)
    ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,save_an_t, term_en)

    """These following arrays aren't the most relevant for this particular 
    type of simulation. Not really needed right now but can tidy up later"""

    np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
    flux_arr.append((ret_flux_frac, pel_pot[p]))
    ener_flux_arr.append((ener_flux, pel_pot[p]))
    lifetime_arr.append((lifetime, pel_pot[p]))

    "Saving analysed data for a singular Bethe calculation"

    np.savetxt(os.path.join(save_an_t, 'terminal_energy_neutral.txt'), term_en)   
    np.savetxt(os.path.join(save_an_t, 'stop_point_neutral.txt'), stop_point, fmt = ('%f'))
    np.savetxt(os.path.join(save_an_t, 'density_neutral.txt'), faux_density)
    np.savetxt(os.path.join(save_an_t, 'real_density_neutral.txt'), real_density)


"Saving compilation of Bethe calculations with varying potentials"
np.savetxt(os.path.join(save_an_t, 'retarded_flux_pot' + str(int(pel_pot[p]))+'.txt'), flux_arr)    
np.savetxt(os.path.join(save_an_t, 'retarded_ener_flux_pot' + str(int(pel_pot[p])) + '.txt'), ener_flux_arr)
np.savetxt(os.path.join(save_an_t, 'pellet_lifetime_retarding_potential_pot' + str(int(pel_pot[p])) +'.txt'), lifetime_arr)
print('success')