import numpy as np 
from gen_var import *
import os
import bethe_analysis_functions as btf 

e_mid, e_bins, MB_norm = dist_calc(e_dist, ener_res, e_bar)
le_mid = len(e_mid)

rp = stop_calc_rp_rc.rp(t)
rc = stop_calc_rp_rc.rc(t) 

min_imp_ener = [] 
particle = 'electron'

arr1 = []
arr2 = []
for t in range(1,lt):
    r = np.arange(0,rc[t] - rp[t],dr) 
    arr1 = btf.stop_term_ener.term_energy(particle,r,t,le)
    arr2, ind = btf.stop_analysis_stop_point.stop_point(arr1,ind,arr2,particle,r)
    arr3 = btf.stop_analysis_particle_density.particle_density(arr2,t,le,ebins,particle)
    min_imp_ener.append(arr1[0,0])
    min_imp_ener.append(t*dt)

arr4 = np.reshape(min_imp_ener, (int(0.5*len(min_imp_ener)),2), 'C')
np.savetxt(os.path.join('Analysed Outputs','GNU_'+particle+'_min_imp_ener_1keV_fs_t.txt'), arr4)
