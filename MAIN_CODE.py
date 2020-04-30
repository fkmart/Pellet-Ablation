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
if style == 'once':
    pot = np.zeros(lr)
    #pot = np.linspace(-2.8*10**-3, 0.0, num = lr, endpoint = 'true')
    savedir_raw = os.path.join(os.getcwd(), 'one_iteration', 'raw_outputs') + os.sep
    savedir_an = os.path.join(os.getcwd(), 'one_iteration', 'analysed_outputs') + os.sep

elif style =='once_charge':
    pot = np.zeros(n_r)
    savedir_raw = os.path.join(os.getcwd(), 'one_iteration_phic', 'raw_outputs') + os.sep
    savedir_an = os.path.join(os.getcwd(), 'one_iteration_phic', 'analysed_outputs') + os.sep

elif style =='many':
    pot = np.zeros(lr)
    savedir_raw = os.path.join(os.getcwd(), 'many_iteration_TEST', 'raw_outputs') + os.sep
    savedir_an = os.path.join(os.getcwd(), 'many_iteration_TEST', 'analysed_outputs') + os.sep

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
acc_elec_dens = np.zeros(rgl) # initial electron density - zero everywhere
pot /= RME*M_fac # normalise electric potential for stopblock calc

if style == 'once':
    p = 0
    fig, ax = plt.subplots()
    for i in range(40, t_end, inc):
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
        #stopblock.stopblock_phi_mod(e_mid, r,i, neut_dens, pot[low:up],save_raw_t) # change this line for new mods
        stopblock.stopblock_phi_mod_rkf(e_mid, r_grid, i, pot, save_raw_t)
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

elif style =='once_charge': 
    save_raw_t = os.path.join(savedir_raw, 't_' + str(i))
    save_an_t = os.path.join(savedir_an, 't_' + str(i))
    if not os.path.exists(save_raw_t):
        os.mkdir(save_raw_t)
    if not os.path.exists(save_an_t):
        os.mkdir(save_an_t)
    low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
    up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
    r_internal = r[low:up]
    mid = int((up-low)/2)
    mid_arr = np.arange(mid,up, 500) 
    p = 10
    neut_dens[low:up] = 0.01*((1.0 + rp[i]**2)/(1.0 + r[low:up]**2))
    for a in range(0, len(sig)):
        print(a)
        for m in mid_arr:
            #Directory Stuff
            phi_max = pel_pot[p]
            path1r = os.path.join(save_raw_t, 'phi_' + str(int(phi_max)))
            path1a = os.path.join(save_an_t, 'phi_' + str(int(phi_max)))
            path2r = os.path.join(path1r, 'sig_{:}'.format(int(100*sig[a])))
            path2a = os.path.join(path1a, 'sig_{:}'.format(int(100*sig[a])))
            path3r = os.path.join(path2r, 'cent_{:}'.format(m))
            path3a = os.path.join(path2a, 'cent_{:}'.format(m))
            if not os.path.exists(path3r):
                if not os.path.exists(path2r):
                    if not os.path.exists(path1r):
                        os.mkdir(path1r)
                    os.mkdir(path2r)
                os.mkdir(path3r)
            if not os.path.exists(path3a):
                if not os.path.exists(path2a):
                    if not os.path.exists(path1a):
                        os.mkdir(path1a)
                    os.mkdir(path2a)
                os.mkdir(path3a)
            #pot = np.linspace(0, pel_pot[p], len(r_internal))
            pot[low:up] = gauss_test_pot.gauss_func(pel_pot[5],sig[a],r_internal[m],r_internal) # using gaussian test function
            pot[:] /= RME*M_fac
            #stopblock.stopblock_phi(e_mid, r,i, neut_dens , pot, path3r)
            stopblock.stopblock_phi_mod_rkf(e_mid, r, i, pot, path3r)
            term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, path3r)
            stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), path3r)
            faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i, len(e_mid), e_bins, particle,path3r,r)
            ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,path3a, term_en)

            #Printing essential information as diagnostic

            print('For peak potential ' + str(pel_pot[p]) + ' the energy flux is ' + str(ener_flux))

            np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
            flux_arr.append((ret_flux_frac, pel_pot[p]))
            ener_flux_arr.append((ener_flux, pel_pot[p]))
            lifetime_arr.append((lifetime, pel_pot[p]))

            "Saving data for a singular Bethe calculation"

            np.savetxt(os.path.join(path3a, 'terminal_energy_pot_test.txt'), term_en)
            np.savetxt(os.path.join(path3a, 'stop_point_pot_test.txt'), stop_point)
            np.savetxt(os.path.join(path3a, 'density_pot_test_mid.txt'), faux_density)
            np.savetxt(os.path.join(path3a, 'real_density_pot_test.txt'), real_density)

elif style =='many':
    i = many_start
    phic = np.zeros(rgl)
    
    while rp[i] > rp_crit and i<80:
        print(i)
        pot = np.zeros(rgl)
        pot[0] = 0.0
        r_grid = np.linspace(rp[i], rc[i], endpoint = 'true', num = rgl)
        low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
        up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
        r_internal = r[low:up]
        phic /= (RME*M_fac)
        stopblock.stopblock_phi_mod_rkf(e_mid,r_grid,i,phic,savedir_raw)

        term_en, ind = baf.stop_analysis_term_ener.term_energy(particle, r, i, le, savedir_raw)
        stop_point = baf.stop_analysis_stop_point.stop_point(term_en,ind, particle, r,i,len(e_mid), savedir_raw)
        faux_density,real_density = baf.stop_analysis_particle_density.particle_density(stop_point,i,len(e_mid),e_bins,particle,savedir_an,r)
        ret_flux_frac, ener_flux, lifetime = baf.stop_analysis_retarded_flux.retarded_flux(i,savedir_an, term_en)

        A = discret_mat.discret(r_grid)
        pot_1 = np.zeros(len(r_grid))
        pot_1[0] = -2.8
        non_dim = r0*r0*e*dens_plas/(epsilon0*1000)
        
        """np.append(flux_arr, (ret_flux_frac, pel_pot[p])) #Append potential dependant arrays with new quantities
        flux_arr.append((ret_flux_frac, pel_pot[p]))
        ener_flux_arr.append((ener_flux, pel_pot[p]))
        lifetime_arr.append((lifetime, pel_pot[p]))"""

        life = life + delta_t

        #shift = j - i # difference between the two "index times" after calculation of new rp
        shift = 10

        #Proposed replacement to commented block beneath
        new_grid = np.linspace(rp[i+shift], rc[i+shift], num = rgl, endpoint = 'true')
        pushed_e_dens = gp.pusher(faux_density, r_grid)
        acc_elec_dens += pushed_e_dens
        phi_before = SOR.SOR(A,pot_1,non_dim*acc_elec_dens,r_grid)
        fig,ax = plt.subplots(figsize = (8,6))
        l1 = ax.plot(r_grid,phi_before, color = 'navy', label = r'$\tilde{\phi}$')
        plt.text(1.0, -0.5, r'$\leftarrow$ to pellet')
        ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
        ax.set_ylabel(r'$\tilde{\phi}$', fontsize = 12, rotation = 0)
        ax2 = ax.twinx()
        ax2.set_ylabel(r'$\tilde{n}$', fontsize =12, rotation = 0)
        l2 = ax2.plot(r_grid,acc_elec_dens, color = 'orange', label = r'$\tilde{n}$')
        lin = l1 + l2 
        lab = [l.get_label() for l in lin]
        plt.legend(lin, lab)
        plt.savefig('phi_before_t' + str(i) + '.png', format = 'png', dpi = 1200)
        plt.show()

        
        e_rem,e_move = elec_transport_push.e_mover(i,shift,r_grid,acc_elec_dens)
        acc_elec_dens = gp.pusher(e_rem,new_grid)
        acc_elec_dens += gp.pusher(e_move, new_grid)
        #Following block is old and comented out until a replacement works
        """
        elec_interp, ind_low, ind_up = com_int.common_interp(r, faux_density[:,0], faux_density[:,1]) # interpolating charge onto all gridpoints
        acc_elec_dens[ind_low:ind_up] = elec_interp[:] + acc_elec_dens[ind_low:ind_up]
        
        "Now move the electrons with the neutrals"
        acc_elec_dens, low_2, up_2 = e_trans.elec_mover(i,r,acc_elec_dens,r[ind_low:ind_up],shift)
        r_domain = r[low_2:up_2] #  shifted spatial range, rp2 to rc2"""
        A = discret_mat.discret(new_grid) #r here needs to be from just in front of the pellet to the cloud at the new time.
        """SOR solver follows - initial guess is zeroes but will be updated to be the old potential with zeroes
        in any extra indices."""
        #phic = SOR.SOR(A, pot[low_2:up_2], acc_elec_dens[low_2:up_2], r_domain)
        non_dim = r0*r0*e*dens_plas/(epsilon0*1000)
        phic = SOR.SOR(A, pot,acc_elec_dens,new_grid)
        fig,ax = plt.subplots(figsize = (8,6))
        l1 = ax.plot(new_grid, phic, label = 'pot')
        plt.text(1.0, -0.5, r'$\leftarrow$ to pellet')
        ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0)
        ax.set_ylabel(r'$\tilde{\phi}$', fontsize = 12)
        ax2 = ax.twinx() 
        ax2.set_ylabel(r'$\tilde{n}$', fontsize = 12, rotation = 0)
        l2 = ax2.plot(new_grid, acc_elec_dens, label = 'dens', color = 'orange')
        lin = l1 + l2 
        lab = [l.get_label() for l in lin]
        plt.legend(lin, lab)
        plt.savefig('phi_after_t' + str(i) + '.png', format = 'png', dpi = 1200)
        plt.show()
        "Saving data for a singular Bethe calculation"
        i+=shift
        np.savetxt(os.path.join(savedir_an, 'terminal_energy_pot_test_t'+str(i)+'.txt'), term_en)   
        np.savetxt(os.path.join(savedir_an, 'stop_point_pot_test_t' + str(i) +'.txt'), stop_point, fmt = ('%f'))
        np.savetxt(os.path.join(savedir_an, 'density_pot_test_t' +str(i) +'.txt'), faux_density)
        np.savetxt(os.path.join(savedir_an, 'real_density_pot_test_t'+str(i) +'.txt'), real_density)
else:
    print("I don't know how you got this far but your style is wrong")
    sys.exit()

"Saving compilation of Bethe calculations with varying potentials"
np.savetxt(os.path.join(save_an_t, 'retarded_flux_pot' + str(int(pel_pot[p]))+'.txt'), flux_arr)    
np.savetxt(os.path.join(save_an_t, 'retarded_ener_flux_pot' + str(int(pel_pot[p])) + '.txt'), ener_flux_arr)
np.savetxt(os.path.join(save_an_t, 'pellet_lifetime_retarding_potential_pot' + str(int(pel_pot[p])) +'.txt'), lifetime_arr)
print('success')