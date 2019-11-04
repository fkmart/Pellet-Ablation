import numpy as np 
import matplotlib.pyplot as plt 
import elec_transport_2 as et 
import iterative_sol as SOR 
import stop_calc_rp_rc as calc 
from gen_var import t, dr, epsilon0
import discret_mat as dm
import datetime as dttm 

rp = calc.calc_rp(t)
rc = calc.calc_rc(t)
shift = 50

print('Spatial resolution is ' + str(dr))
load_dir = '/home/kyle/Documents/Python_Code/stop_code/one_iteration/analysed_outputs'

t_ind = 50
elec_dens_1 = np.loadtxt(load_dir + '/real_density_t'+str(t_ind)+'.txt')
r = np.arange(0, rc[-1], dr)

elec_dens_2, r_plot_2, ed_1, r_plot_1 = et.elec_mover(t_ind, r, elec_dens_1[:,0], elec_dens_1[:,1], shift)

#Now solve Poisson for these densities 
rp1 = rp[t_ind]

rp2 = rp[t_ind+shift]
rc2 = rc[t_ind + shift]

indl1 = next(p[0] for p in enumerate(r) if p[1] > rp1)

indl2 = next(p[0] for p in enumerate(r) if p[1] > rp2)
indu2 = next(p[0] for p in enumerate(r) if p[1] > rc2)

r_sor = r[indl2:indu2]

A = dm.discret(r_sor)
phi = np.zeros(len(r_sor))
phi[0] = -3

time1 = dttm.datetime.now()

K = 1.0/(epsilon0*1000)

x = SOR.SOR(A, phi, -K*elec_dens_2,r)

time2 = dttm.datetime.now()


#Check the profile has actually shifted 

rmax1 = r_plot_1[np.argmax(ed_1)]
rmax2 = r_plot_2[np.argmax(elec_dens_2)]


#print('Time for solver at t = ' + str(t_ind + shift) + ' is : ' + str(time2-time1))
fig, ax = plt.subplots()
ax.plot(elec_dens_1[:,1]+r[indl1 -1], elec_dens_1[:,0])
ax.plot(r_plot_1,ed_1 )
ax.plot(r_plot_2,elec_dens_2)
#ax.axvline(rmax1, linestyle = '--')
#ax.axvline(rmax2, linestyle = '--')
ax2 = ax.twinx()
ax2.plot(r_sor, x, color = 'firebrick', linestyle = '--')
#ax.set_yscale('log')
plt.grid(axis = 'both')
plt.show()

print('Did it work?')