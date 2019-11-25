"""This code will use an RK4 solver to solve the CSDA equation from Sletzer and Berger to determine the distribution
of electrons across the cloud due to inelastic collisions with neutral diatomic hydrogen. The distribution of charge
is then used to determine the self-field according to Poisson's equation by using an SOR method. The incident 
energy on the pellet surface then determines the """
import numpy as np
from gen_var import rp, rc, t_start,lt, lr, dr, lp, pel_pot, cloud_pot, style , many_start, t_end, inc, p_inc, rp_crit, delta_t, lp , sig, r, n_r
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

neut_dens = np.zeros(n_r) # neutral density array of equal length, only some values will be non-zero

lr = len(r)

#Important to print the style being used here so you don't totally miss it
print('The scenario tested here is: ' + str(style))

"Must now establish how the simulation will proceed with choice of style"
if style == 'once':
    pot = np.zeros(lr)
    savedir_raw = os.path.join(os.getcwd(), 'one_iteration', 'raw_outputs') + os.sep
    savedir_an = os.path.join(os.getcwd(), 'one_iteration', 'analysed_outputs') + os.sep
    #savedir_raw = os.getcwd() + "/one_iteration/raw_outputs/"
    #savedir_an = os.getcwd() + "/one_iteration/analysed_outputs/"
elif style =='once_charge':
    pot = np.zeros(n_r)
    savedir_raw = os.path.join(os.getcwd(), 'one_iteration_phic', 'raw_outputs') + os.sep
    savedir_an = os.path.join(os.getcwd(), 'one_iteration_phic', 'analysed_outputs') + os.sep
    #savedir_raw = os.getcwd() + "/one_iteration_phic/raw_outputs/"
    #savedir_an = os.getcwd() + "/one_iteration_phic/analysed_outputs/"
elif style =='many':
    pot = np.zeros(lr)
    savedir_raw = os.path.join(os.getcwd(), 'many_iteration', 'raw_outputs') + os.sep
    savedir_an = os.path.join(os.getcwd(), 'many_iteration', 'analysed_outputs') + os.sep
    #savedir_raw = os.getcwd() + "/many_iteration/raw_outputs/"
    #savedir_an = os.getcwd() + "/many_iteration/analysed_outputs/"
else:
    print ('Must select an appropriate simulation style - see gen_var for details')
    sys.exit()

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
acc_elec_dens = np.zeros(n_r) # initial electron density - zero everywhere
pot /= RME*M_fac # normalise electric potential for stopblock calc

if style == 'once':
    p = 0
    fig, ax = plt.subplots()
    for i in range(t_start, t_end, inc):
        print(i)
        low = next(p[0] for p in enumerate(r) if p[1] > rp[i]) # finds first index with r > r_p
        up = next(p[0] for p in enumerate(r) if p[1] > rc[i]) # finds first index with r > r_c
        r_internal = r[low:up] 
        neut_dens[low:up] = 0.01*((1.0 + rp[i]**2)/(1.0 + r[low:up]**2))
        pot = np.zeros(n_r)
        stopblock.stopblock_phi(e_mid, r,i, neut_dens, pot,savedir_raw)

        term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, savedir_raw)
        stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), savedir_raw)
        faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,p,r_internal[0],r)
        ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,savedir_an, term_en)

        """These following arrays aren't the msot relevant for this particular 
        type of simulation. Not really needed right now but can tidy up later"""

        np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
        flux_arr.append((ret_flux_frac, pel_pot[p]))
        ener_flux_arr.append((ener_flux, pel_pot[p]))
        lifetime_arr.append((lifetime, pel_pot[p]))

        "Saving analysed data for a singular Bethe calculation"

        np.savetxt(os.path.join(savedir_an, 'terminal_energy_pot_test_t'+str(i)+'.txt'), term_en)   
        np.savetxt(os.path.join(savedir_an, 'stop_point_pot_test_t' + str(i)+'.txt'), stop_point, fmt = ('%f'))
        np.savetxt(os.path.join(savedir_an, 'density_pot_test_t' +str(i) +'.txt'), faux_density)
        np.savetxt(os.path.join(savedir_an, 'real_density_pot_test_t'+str(i) +'.txt'), real_density)

elif style =='once_charge': 
    low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
    up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
    r_internal = r[low:up]
    mid = int((up-low)/2)
    neut_dens[low:up] = 0.01*((1.0 + rp[i]**2)/(1.0 + r[low:up]**2))
    for a in range(0, len(sig)):
        print(a)
        for p in range(0, lp, p_inc):
            #pot = np.linspace(0, pel_pot[p], len(r_internal))
            pot[low:up] = gauss_test_pot.gauss_func(pel_pot[p],sig[a],r_internal[mid],r_internal) # using gaussian test function
            pot[:] /= RME*M_fac
            stopblock.stopblock_phi(e_mid, r,i, neut_dens , pot, savedir_raw)

            term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, savedir_raw)
            stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), savedir_raw)
            faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,p,r_internal[0],r)
            ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,savedir_an, term_en)

            #Printing essential information as diagnostic

            print('For peak potential ' + str(pel_pot[p]) + ' the energy flux is ' + str(ener_flux))

            np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
            flux_arr.append((ret_flux_frac, pel_pot[p]))
            ener_flux_arr.append((ener_flux, pel_pot[p]))
            lifetime_arr.append((lifetime, pel_pot[p]))

            "Saving data for a singular Bethe calculation"

            np.savetxt(os.path.join(savedir_an, 'terminal_energy_pot_test_t'+str(i)+'pot'+str(pel_pot[p])+'sig'+ str(sig[a]) +'.txt'), term_en)   
            np.savetxt(os.path.join(savedir_an, 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'sig' + str(sig[a])+'.txt'), stop_point, fmt = ('%f'))
            np.savetxt(os.path.join(savedir_an, 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'sig' + str(sig[a]) +'.txt'), faux_density)
            np.savetxt(os.path.join(savedir_an, 'real_density_pot_test_t'+str(i) +'pot'+str(pel_pot[p])+'sig' + str(sig[a])+'.txt'), real_density)

elif style =='many':
    i = many_start
    while rp[i] > rp_crit:
        
        low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
        up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
        r_internal = r[low:up]
        stopblock.stopblock_phi(e_mid, r,i, neut_dens, pot, savedir_raw)

        term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, savedir_raw)
        stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), savedir_raw)
        faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,savedir_raw, r_internal[0],r)
        ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,savedir_an, term_en)

        """np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
        flux_arr.append((ret_flux_frac, pel_pot[p]))
        ener_flux_arr.append((ener_flux, pel_pot[p]))
        lifetime_arr.append((lifetime, pel_pot[p]))"""

        life = life + delta_t

        #shift = j - i # difference between the two "index times" after calculation of new rp
        shift = 1

        elec_interp, ind_low, ind_up = com_int.common_interp(r, real_density[:,0], real_density[:,1]) # interpolating charge onto all gridpoints
        acc_elec_dens[ind_low:ind_up] = elec_interp[:] + acc_elec_dens[ind_low:ind_up]
        
        "Now move the electrons with the neutrals"
        acc_elec_dens, low_2, up_2 = e_trans.elec_mover(i,r,acc_elec_dens,r[ind_low:ind_up],shift)
        r_domain = r[low_2:up_2] #  shifted spatial range, rp2 to rc2
        A = discret_mat.discret(r_domain) #r here needs to be from just in front of the pellet to the cloud at the new time.
        """SOR solver follows - initial guess is zeroes but will be updated to be the old potential with zeroes
        in any extra indices."""
        #phic = SOR.SOR(A, pot[low_2:up_2], acc_elec_dens[low_2:up_2], r_domain)

        "Saving data for a singular Bethe calculation"

        np.savetxt(os.path.join(savedir_an, 'terminal_energy_pot_test_t'+str(i)+'.txt'), term_en)   
        np.savetxt(os.path.join(savedir_an, 'stop_point_pot_test_t' + str(i) +'.txt'), stop_point, fmt = ('%f'))
        np.savetxt(os.path.join(savedir_an, 'density_pot_test_t' +str(i) +'.txt'), faux_density)
        np.savetxt(os.path.join(savedir_an, 'real_density_pot_test_t'+str(i) +'.txt'), real_density)

else:
    print("I don't know how you got this far but your style is wrong")
    sys.exit()

"Saving compilation of Bethe calculations with varying potentials"
np.savetxt(os.path.join(savedir_an, 'retarded_flux_pot_test_t'+str(i)+'pot' + str(pel_pot[p])+'.txt'), flux_arr)    
np.savetxt(os.path.join(savedir_an, 'retarded_ener_flux_pot_test'+str(i)+ 'pot' + str(pel_pot[p]) + '.txt'), ener_flux_arr)
np.savetxt(os.path.join(savedir_an, 'pellet_lifetime_retarding_potential_pot_test' + str(i)+ 'pot' + str(pel_pot[p]) +'.txt'), lifetime_arr)
print('success')