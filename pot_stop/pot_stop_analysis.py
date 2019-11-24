import numpy as np
from electron import e_bar, e_dist, ener_res, dist_calc 
from gen_var import rp, rc , t_start, sig, pel_pot , lp, I 
import gauss_test_pot as gtp
import bethe_analysis_functions as baf 
import os 

direc = os.getcwd() + '/pot_stop'
subdir_load = '/data/'
subdir_save = '/analysed_data/'

e_mid, e_bins, mb = dist_calc(e_dist, ener_res,e_bar)

r_int = np.linspace(rp[t_start], rc[t_start], 602)

for a in range(0, len(sig)):
    for p in range(1, lp):
        term_ener, indices = baf.stop_analysis_term_ener.term_energy(r_int,p,a, len(e_mid), direc+subdir_load)
        stop_point = baf.stop_analysis_stop_point.stop_point(indices, r_int, p,a, len(e_mid), direc + subdir_load)
        faux_dens, real_dens = baf.stop_analysis_particle_density.particle_density(stop_point,t_start, len(e_mid),e_bins, r_int)
        y = format(sig[a], '.2f')
        np.savetxt(direc + subdir_save + 'term_ener_peak_' + str(int(pel_pot[p])) + 'sig_' + y + '.txt', term_ener)
        np.savetxt(direc + subdir_save + 'stop_point_peak' + str(int(pel_pot[p])) + 'sig_' + y + '.txt', stop_point)
        np.savetxt(direc + subdir_save + 'elec_bins_peak' + str(int(pel_pot[p])) + 'sig_' + y + '.txt', faux_dens)
        np.savetxt(direc + subdir_save + 'elec_dens_peak' + str(int(pel_pot[p])) + 'sig_' + y + '.txt', real_dens)
print('I think you might have done it')