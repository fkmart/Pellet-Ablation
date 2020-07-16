"""This code will use an RK4 solver to solve the CSDA equation from Sletzer and Berger to determine the distribution
of electrons across the cloud due to inelastic collisions with neutral diatomic hydrogen. The distribution of charge
is then used to determine the self-field according to Poisson's equation by using an SOR method. The incident 
energy on the pellet surface then determines the """
import numpy as np
from gen_var import rp_hr, rc_hr, t_start,lt, lr, dr, lp, pel_pot, cloud_pot, style,r0, dens_plas,e,epsilon0
from gen_var import many_start, t_end, inc, p_inc, rp_crit, delta_t, sig, rgl, r_grid
from gen_var import N_0_sca, pel_dens_numb, phi_plas,dt_hr, t_low_hr,t_upper_hr,tf, phi_p, sheath_pot, sheath_x
import os 
import stopblock #function to bethe stop electrons
import MB_calc # function to determine EEDF
from electron import e_dist, ener_res, e_bar, KE_top, KE_bot, RME, M_fac #importing essential variables and functions
import electron as part
import bethe_analysis_functions as baf # funcitons to analyse the data
import sys #for exiting the code and printing an "error message"
import discret_mat # makes discretisation matrix for SOR
import iterative_sol as SOR #SOR solver
#import elec_transport_2 as e_trans #transport module to move electrons based on relative change in H2 densities
#import common_interpolator as com_int #interpolater to move from one grid to another
import matplotlib.pyplot as plt # plotting module - not to be used in these simulations but present for testing
import grid_pusher as gp
#import elec_transport_push 
from TRANSPORT_FUNC_FINAL import transport
import scipy.interpolate as spint
import ablation
import find_nearest as fn
import norming_density as norm
import remove_jumps_func as rj
import number_zero

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

push_e_dens = np.zeros(rgl)
#lr = len(r)

#Important to print the style being used here so you don't totally miss it
print('The scenario tested here is: ' + str(style))

"Must now establish how the simulation will proceed with choice of style"

pot = np.zeros(rgl)


#######################################################################
TYPE = 'sheath'
#######################################################################


raw_data = os.path.join(os.getcwd(), 'many_iteration',TYPE, str(int(e_bar)) +'eV', 'raw_outputs') + os.sep
an_data = os.path.join(os.getcwd(), 'many_iteration',TYPE,str(int(e_bar))+'eV', 'analysed_outputs') + os.sep

"SAVE INITIAL PARAMETERS"
params = ['t_low_hr (nd)','t_upper_hr (nd)','dt_hr (nd)','delta_t (nd)','tf (s)','bin_width (eV)','e_bar (eV)']
vals = [str(t_low_hr),str(t_upper_hr),str(dt_hr),str(delta_t),str(tf),str(ener_res), str(e_bar)]
param_save = [params,vals]
np.savetxt(an_data +'params.txt',param_save, fmt = '%s', delimiter = ',')

"Establish the new arrays"
flux_arr = []
ener_flux_arr = []
lifetime_arr = []
acc_elec_dens = np.zeros(rgl) # initial electron density - zero everywhere
phic = np.zeros(rgl)

i = many_start
if TYPE == 'poisson':
    phic[0] = 0.0
else:
    phic[0] = 0.0


N_p_n = ((4.0*np.pi/3.0)*pel_dens_numb*r0**3 * rp_hr[i]**3)/N_0_sca

count = 0
shift_save = []
fractions = []
time_indices = []
rp_rc_arr = []
phi_elec = np.zeros(rgl)
#r_grid = np.linspace(rp_hr[i], rc_hr[i],num = rgl, endpoint = 'true')
ind_low_r = fn.find_nearest(r_grid,rp_hr[i])
ind_up_r = fn.find_nearest(r_grid,rc_hr[i])

if TYPE =='sheath':
    sheath_x = np.abs(sheath_x)/r0 + rp_hr[i]
    phi_sheath = spint.pchip_interpolate(sheath_x, sheath_pot, r_grid)
    phi_sheath[:ind_low_r] = 0.0
    for a in range(0, len(phi_sheath)):
        if phi_sheath[a] > 0.0:
            phi_sheath[a] = 0.0
    phic = np.copy(phi_sheath)/(RME*M_fac)


N_prev = 0.0
rp_rc_arr.append((rp_hr[i],rc_hr[i]))
while rp_hr[i] > rp_crit and count < 20:
    print('Iteration number is ' + str(count))
    print(i)    
    time_indices.append(i)
    if TYPE =='neutral':
        phic = np.zeros(rgl)
    else:
        pass
    path_raw_new = os.path.join(raw_data, str(count))
    path_an_new = os.path.join(an_data,str(count))
    if os.path.exists(path_raw_new) and os.path.exists(path_an_new):
        pass
    else:
        os.mkdir(path_raw_new)
        os.mkdir(path_an_new)
    savedir_raw = path_raw_new
    savedir_an = path_an_new
    #r_grid = np.linspace(rp_hr[i], rc_hr[i], endpoint = 'true', num = rgl)
    #low = next(p[0] for p in enumerate(r) if p[1] > rp_hr[i])
    #up = next(p[0] for p in enumerate(r) if p[1] > rc_hr[i])
    #r_internal = r[low:up]
    stopblock.stopblock_phi_mod_rkf(e_mid,r_grid,i,phic,savedir_raw)

    term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, i, le, savedir_raw)
    stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle,i,len(e_mid), savedir_raw)
    faux_density = baf.stop_analysis_particle_density.particle_density(stop_point,i,len(e_mid),e_bins,particle,savedir_an)
    norm_density = norm.norm(i,faux_density, r_grid[ind_low_r:ind_up_r])
    #ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i, term_en)
    N_loss_n, ener_flux, lifetime, fraction = ablation.fluxes(i,faux_density, term_en)
    N_loss_n += N_prev
    rp_new = ablation.new_rp(N_loss_n, N_p_n,i)
    N_p_n = rp_new**3
    life += delta_t

    #shift = j - i # difference between the two "index times" after calculation of new r
    shift = fn.find_nearest(rp_hr,rp_new) - i
    if shift ==0:
        N_prev = np.copy(N_loss_n)
    else:
        N_prev = 0.0
#    A = discret_mat.discret(r_grid[ind_low:ind_up])
    #np.append(flux_arr, (ret_flux_frac, i)) #Append potential dependant arrays with new quantities
    #flux_arr.append((ret_flux_frac, i))
    ener_flux_arr.append((ener_flux, i))
    lifetime_arr.append((lifetime, i))
    fractions.append((fraction, i))
    #Proposed replacement to commented block beneath
    """new_grid = np.linspace(rp[i+shift], rc[i+shift], num = rgl, endpoint = 'true')
    pushed_e_dens = gp.pusher(faux_density, r_grid)
    acc_elec_dens += pushed_e_dens
    phi_before = SOR.SOR(A,pot_1,non_dim*acc_elec_dens,r_grid)"""
    #e_rem,e_move = elec_transport_push.e_mover(i,shift,r_grid,acc_elec_dens)
    #new_grid = np.linspace(rp_hr[i+shift], rc_hr[i+shift], num = rgl, endpoint = 'true')
    
    #g = spint.interp1d(faux_density[:,1],faux_density[:,0], kind = 'cubic', fill_value = 'extrapolate' )
    #push_interp_elec_dens = g(r_grid)
    #push_elec_dens = gp.pusher(faux_density, r_grid) # interpolate points here to the grid
    #acc_elec_dens += push_interp_elec_dens
    acc_elec_dens[ind_low_r:ind_up_r] += np.flip(norm_density,axis = 0)
    e_rem, e_trans = transport(i,shift, r_grid, acc_elec_dens) # transported fractions of density defined
    acc_elec_dens = gp.pusher(e_rem,r_grid) # remaining fraction pushed to new grid
    ind_up, ind_low = number_zero.index_crit(acc_elec_dens) #POSSIBLY REMOVE THIS
    acc_elec_dens += gp.pusher(e_trans, r_grid) #transported fraction pushed to new grid
    "Now need to remove any jumps due to erroneous pushing"
    #acc_elec_dens = rj.smoother(r_grid, acc_elec_dens, ind_up, ind_low)#POSSIBLY REMOVE THIS
    #acc_elec_dens *= (count +1)/np.sum(acc_elec_dens)
    
    ind_low_r = fn.find_nearest(r_grid,rp_hr[i+shift])
    ind_up_r = fn.find_nearest(r_grid, rc_hr[i+shift])
    #Following block is old and comented out until a replacement works
    """
    elec_interp, ind_low, ind_up = com_int.common_interp(r, faux_density[:,0], faux_density[:,1]) # interpolating charge onto all gridpoints
    acc_elec_dens[ind_low:ind_up] = elec_interp[:] + acc_elec_dens[ind_low:ind_up]
        
    "Now move the electrons with the neutrals"
    acc_elec_dens, low_2, up_2 = e_trans.elec_mover(i,r,acc_elec_dens,r[ind_low:ind_up],shift)
    r_domain = r[low_2:up_2] #  shifted spatial range, rp2 to rc2"""
    A = discret_mat.discret(r_grid[ind_low_r:ind_up_r]) #r here needs to be from just in front of the pellet to the cloud at the new time.
    """SOR solver follows - initial guess is zeroes but will be updated to be the old potential with zeroes
    in any extra indices."""
    #phic = SOR.SOR(A, pot[low_2:up_2], acc_elec_dens[low_2:up_2], r_domain)
    #if count != 0:
    #    phic *= RME*M_fac/phi_plas
    #else:
    #    pass
    non_dim = r0*r0*e*dens_plas/(epsilon0*phi_plas)
    phi_elec[ind_low_r:ind_up_r] = SOR.SOR(A, phi_elec[ind_low_r:ind_up_r],non_dim*acc_elec_dens[ind_low_r:ind_up_r],r_grid[ind_low_r:ind_up_r])
    if TYPE =='sheath':
        x_thing = sheath_x - rp_hr[i] + rp_hr[i + shift]
        phi_sheath = spint.pchip_interpolate(sheath_x - rp_hr[i] + rp_hr[i + shift],sheath_pot, r_grid)/phi_plas
        phi_sheath[:ind_low_r] = 0.0
        for a in range(0, len(phi_sheath)):
            if phi_sheath[a] > 0.0:
                phi_sheath[a] = 0.0
        phic = np.copy(phi_elec) 
        phic += phi_sheath
    else:
        phic = np.copy(phi_elec) 
    phic[:ind_low_r]= 0.0
    phic[ind_up_r:] = 0.0
    "Saving data for a singular Bethe calculation"
    
    np.savetxt(os.path.join(savedir_an, 'terminal_energy_t'+str(i)+'.txt'), term_en)   
    np.savetxt(os.path.join(savedir_an, 'stop_point_t' + str(i) +'.txt'), stop_point, fmt = ('%f'))
    np.savetxt(os.path.join(savedir_an, 'density_t' +str(i) +'.txt'), faux_density)
    np.savetxt(os.path.join(savedir_an, 'real_density_t'+str(i) +'.txt'), norm_density)
    np.savetxt(os.path.join(savedir_an, 'potential_t' + str(i) + '.txt'), np.asarray([phic*phi_plas,r_grid]))
    np.savetxt(os.path.join(savedir_an, 'accummulated_density_t' + str(i) + '.txt'), np.asarray([acc_elec_dens,r_grid]))
    i+=shift
    count += 1
    shift_save = np.append(shift_save,shift)
    rp_rc_arr.append((rp_hr[i],rc_hr[i]))
    phic *= phi_plas/(RME*M_fac)
"Saving compilation of Bethe calculations with varying potentials"
#np.savetxt(os.path.join(an_data, 'retarded_flux.txt'), flux_arr)    
np.savetxt(os.path.join(an_data,'time_indices.txt'),time_indices)
np.savetxt(os.path.join(an_data, 'retarded_ener_flux.txt'), ener_flux_arr)
np.savetxt(os.path.join(an_data, 'pellet_lifetime' +'.txt'), lifetime_arr)
np.savetxt(os.path.join(an_data, 'time_shifts.txt'), shift_save)
np.savetxt(os.path.join(an_data, 'impacting_fractions.txt'), fractions)
np.savetxt(os.path.join(an_data, 'rp_rc.txt'), rp_rc_arr)
print('success')