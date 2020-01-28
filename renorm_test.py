import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import lp, p_inc, pel_pot, rp, rc , t_start, sig , r, n_r
import TEST_stop_analysis_norming_elec_density as stned 
import electron 
from electron import e_bar, e_dist, ener_res 
import romberg as ro
import FUNC_renorm_dens as fnd

e1,e_bins, e2 = electron.dist_calc(e_dist, ener_res, e_bar)
direc = os.getcwd() 
load_dir = os.path.join(direc, 'one_iteration_phic', 'analysed_outputs') + os.sep
#load_dir = direc + '/one_iteration_phic/analysed_outputs/'
savedir = os.path.join(direc,'pictures') + os.sep
i = t_start 

k = sig[-1]
rp = rp[i]
fig, ax = plt.subplots(figsize = (10,8))

#for a in range(0, len(sig)):
#    k = sig[a]
integ_arr = []
"""for p in range(0, lp-1, 1):
    file = np.loadtxt(load_dir + 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'sig' + str(k) +'.txt')
    arr2 = np.loadtxt(load_dir + 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'sig' + str(k)+'.txt')
    real_dens, integrated,rdf,fd_int = stned.renorm_dens(file[:,0],file[:,1],e_bins,arr2,i,r )
    integ_arr = np.append(integ_arr, (integrated, pel_pot[p]))
    #ax.plot(real_dens[:,1], real_dens[:,0], label = r'$\phi = $' + str(pel_pot[p]) + 'V')
    #ax.plot(rdf[:750,1],rdf[:750,0], label = r'$\phi = $' + str(pel_pot[p]))
    ax.plot(file[:,0], file[:,1], label = r'$\phi = $' + str(pel_pot[p]) + 'V') 
    #ax.plot(fd_int[:,1], fd_int[:,0], label = r'$\phi = $' + str(pel_pot[p]))
    #print(ro.romberg_samp(fd_int[:,0], fd_int[:,1]))
    "To test the difference between interpolated bins and raw bin data"
    s = sum(file[:,1])

    #fd_int[:,0] = fd_int[:,0]/(sum(fd_int[:,0]))
    #fd_int[:,0] *= s 
    #ax.plot(file[:,0], file[:,1], label = 'raw - ' + str(pel_pot[p]))
    #ax.plot(fd_int[:,1], fd_int[:,0], label = 'interp - ' + str(pel_pot[p]))
    print(str(integrated) + ' and the peak potential is ' + str(pel_pot[p]))"""
for p in range(0, lp-1, 1):
    file = np.loadtxt(load_dir + 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'sig' + str(k) +'.txt')
    arr2 = np.loadtxt(load_dir + 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'sig' + str(k)+'.txt')
    bins_full, dens_full,r_full = fnd.renorm(file[:,0], file[:,1])
    ax.plot(r_full, dens_full, label = r'$\phi_t = $' + str(pel_pot[p]) + 'V')

#integ_arr = np.reshape(integ_arr, (int(len(integ_arr)*0.5),2), 'C')
#ax.plot(integ_arr[:,1], integ_arr[:,0], label = r'$\sigma = $' + str(k))
plt.legend(ncol = 3, loc = 'lower right')
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel(r'$\mathrm{Density}/\mathrm{m}^{-3}$')
ax.set_yscale('log')
#ax_2 = ax.twinx().twiny()
#ax_2.scatter(integ_arr[:,1], integ_arr[:,0], color = 'k', marker = 'x' )
#plt.savefig(load_dir + 'integrated_val_gauss_pot_sig_all.png', format = 'png', dpi = 1200)
"Following code is to compare interpolation to raw data"
#p = 10
#file = np.loadtxt(load_dir + 'density_pot_test_t' +str(i) +'pot'+str(pel_pot[p])+'sig' + str(k) +'.txt')
#arr2 = np.loadtxt(load_dir + 'stop_point_pot_test_t' + str(i) +'pot'+str(pel_pot[p])+'sig' + str(k)+'.txt')
#real_dens, integrated,rdf = stned.renorm_dens(file[:,0],file[:,1],e_bins,arr2,i,r )
#ax.plot(real_dens[:,1], real_dens[:,0], label='raw')
#ax.plot(rdf[:,1], rdf[:,0], label = 'interpolated')
#plt.legend()
plt.savefig(savedir + 'dens_fix.png', format = 'png', dpi = 1200)
plt.show()