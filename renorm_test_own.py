import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import lp, p_inc, pel_pot, rp, rc , t_start, sig , r, n_r
import TEST_stop_analysis_norming_elec_density as stned 
import electron 
from electron import e_bar, e_dist, ener_res 

e1,e_bins, e2 = electron.dist_calc(e_dist, ener_res, e_bar)
direc = os.getcwd() 
load_dir = os.path.join(direc, 'one_iteration_phic', 'analysed_outputs') + os.sep
#load_dir = direc + '/one_iteration_phic/analysed_outputs/'
i = t_start 

k = sig[-2]
rp = rp[i]
fig, ax = plt.subplots()
integ_arr = [] 
integ_arr_2 = []
for p in range(0, lp, p_inc):
    file = np.loadtxt(load_dir + 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'sig' + str(k) +'.txt')
    arr2 = np.loadtxt(load_dir + 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'sig' + str(k)+'.txt')
    real_dens, integrated,rdf,fd_int, integ_2 = stned.renorm_dens(file[:,0],file[:,1],e_bins,arr2,i,r )
    integ_arr = np.append(integ_arr, (integrated, pel_pot[p]))
    integ_arr_2 = np.append(integ_arr_2, (integ_2, pel_pot[p]))
integ_arr = np.reshape(integ_arr, (int(len(integ_arr)*0.5),2), 'C')
integ_arr_2 = np.reshape(integ_arr_2, (int(len(integ_arr_2)*0.5),2), 'C')
ax.plot(integ_arr[:,1], integ_arr[:,0], label =  'scipy')
ax.scatter(integ_arr[:,1], integ_arr[:,0], label = 'homemade', marker = 'x')

plt.legend()
ax.set_xlabel('Peak Retarding Potential/V')
ax.set_ylabel('Line Integrated Value')
#ax.set_yscale('log')
plt.show()