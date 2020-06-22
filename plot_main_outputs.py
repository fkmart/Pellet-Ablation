import numpy as np
import matplotlib.pyplot as plt 
import os 
from gen_var import lp, pel_pot, p_inc, rp, rc 
import gauss_test_pot

#load_dir = os.getcwd() + '/one_iteration_phic/analysed_outputs/'
load_dir = os.path.join(os.getcwd(), 'one_iteration', 'analysed_outputs') + os.sep
#save_dir = '/home/kyle/Documents/thesis/Local_figures/'
save_dir = os.path.join(os.getcwd(), 'pictures') + os.sep
i = 50

x_res = 0.01 
n = 1
n_r = 2.0**n
while (rc[-1]/(2.0**(n)) > x_res) : 
    n += 1

n_r = 2**n + 1

r = np.linspace(0, rc[-1], n_r) # romberg grid defined

low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
up = next(p[0] for p in enumerate(r) if p[1] > rc[i])
r_internal = r[low:up]
mid = int(up/2)
k = 2.0
pot = gauss_test_pot.gauss_func(1000,k,r_internal[mid],r_internal) # using gaussian test function
print('Potential at cloud surface is ' + str(pot[-1]))
print('Potential at pellet surface is ' + str(pot[0]))
fig, ax = plt.subplots()

plot_arr = [] 
for p in range(0, lp, p_inc +1):
    file = np.loadtxt(load_dir + 'density_pot_test_t'+str(i) +'pot'+str(pel_pot[p])+'.txt')
    #file_min = np.loadtxt(load_dir + 'terminal_energy_pot_test_t50pot'+str(pel_pot[p]) + '.txt')
    #ax.plot(file[:,0], file[:,1], label = r'$\phi_{\mathrm{max}} = ' + str(pel_pot[p]) + 'V $')
    #plot_arr.append(file_min[0,0])
    #plot_arr.append(pel_pot[p])
#plot_arr = np.asarray(plot_arr) 
#plot_arr = np.reshape(plot_arr, (int(len(plot_arr)/2),2),'C')
#file = np.loadtxt(load_dir + 'retarded_flux_pot_test_t'+str(i)+'pot-2000.0.txt')

file1 = np.loadtxt(load_dir + 'terminal_energy_pot_test_t50pot-0.0.txt')
file2 = np.loadtxt(load_dir + 'terminal_energy_pot_test_t50pot-2000.0.txt')
init_ener = file1[-100:,0]
imp1 = file1[-100:,1]
imp2 = file2[-100:,1]
imp_diff = imp1 - imp2

ax.plot(init_ener, imp_diff)
ax.set_ylabel('Impacting Energy Difference/eV')
ax.set_xlabel('Initial Electron Energy/eV')
#plt.legend()
#ax2 = ax.twinx()
#ax2.plot(r_internal,pot, linestyle = '--') 
plt.savefig(save_dir + 'gauss_prof_k' + str(k) +'_pot2000_impenerdiff.png', format = 'png', dpi = 1200)
plt.show()