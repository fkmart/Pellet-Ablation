import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import lp, p_inc, pel_pot, rp, rc , t_start, sig , r, n_r
import stop_analysis_norming_elec_density as stned 
import electron 
from electron import e_bar, e_dist, ener_res 

e1,e_bins, e2 = electron.dist_calc(e_dist, ener_res, e_bar)
direc = os.getcwd() 
load_dir = direc + '/one_iteration_phic/analysed_outputs/'
i = t_start 

k = sig[-3]
rp = rp[i]

fig, ax = plt.subplots()
for p in range(0, lp, p_inc):
    file = np.loadtxt(load_dir + 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'sig' + str(k) +'.txt')
    arr2 = np.loadtxt(load_dir + 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'sig' + str(k)+'.txt')
    real_dens = stned.renorm_dens(file[:,0],file[:,1],e_bins,arr2,i,r )
    #integ_arr = np.append(integ_arr, (integrated, pel_pot[p]))
    #ax.plot(real_dens[:,1], real_dens[:,0], label = r'$\phi = $' + str(pel_pot[p]))
    #ax.plot(r[:750],rdf[:750], label = r'$\phi = $' + str(pel_pot[p]))
    #ax.plot(file[:,0], file[:,1], label = r'$\phi = $' + str(pel_pot[p])) 
#integ_arr = np.reshape(integ_arr, (int(len(integ_arr)*0.5),2), 'C')
#ax.plot(integ_arr[:,1], integ_arr[:,0], label = r'$\sigma = $' + str(k))
#ax.set_xlabel('Peak Retarding Potential/V')
#ax.set_ylabel('Line Integrated Value')
#ax.set_yscale('log')
#plt.legend()
#plt.savefig(load_dir + 'integrated_val_gauss_pot_sig_all.png', format = 'png', dpi = 1200)
plt.show()
