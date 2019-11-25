import numpy as np
import matplotlib.pyplot as plt
from gen_var import dr,dt , t, eps 
import stop_calc_rp_rc
import scipy.integrate as spint
import scipy.interpolate as spit
import os

mydir = os.path.join(os.getcwd(), 'Local_figures') + os.sep
rp = stop_calc_rp_rc.calc_rp(t)
rc = stop_calc_rp_rc.calc_rc(t)
rcd = stop_calc_rp_rc.rdot(t,rc,rp)
 
t1 = 50
t2 =150
delta_t = t2 - t1

#common grid
dr = 0.001
r = np.arange(0,rc[-1], dr)

"for t1"
rp1 = rp[t1] #pellet radius
rc1 = rc[t1] #cloud radius
rcd1 = rcd[t1] #cloud radius expansion speed
#r1 = np.arange(rc1,rp1, -dr) 
dens1 = eps*((1.0 + rp1**2)/(1.0 + r**2))

"for t2" # same as above lines
rp2 = rp[t2] 
rc2 = rc[t2]
rcd2 = rcd[t2] 
#r2 = np.arange(rc2, rp2, -dr)
dens2 = eps*((1.0 + rp2**2)/(1.0 + r**2))

x = r - rp1
#following four lines calculate the index corresponding to rp1,rp2,rc1 and rc2
ind_low_1 = (np.abs(r - rp1)).argmin() +1
ind_up_2 = (np.abs(r - rc2)).argmin() +1
ind_up_1 = (np.abs(r - rc1)).argmin()+1
ind_low_2 = (np.abs(r - rp2)).argmin() + 1

#zero densities outside determined limits
dens1[:ind_low_1] = 0.0
dens2[:ind_low_2] = 0.0 
dens1[ind_up_1:] = 0.0
dens2[ind_up_2:] = 0.0 

net_change = (dens1[ind_low_1:ind_up_1] - dens2[ind_low_1:ind_up_1]) #calculates netchange in density for common points

#plotting stuff

r_plot1 = r[ind_low_1:ind_up_1]
r_plot2 = r[ind_low_2:ind_up_2]
dens1_plot = dens1[ind_low_1:ind_up_1]
dens2_plot = dens2[ind_low_2:ind_up_2]
new_dens = dens1[ind_low_1:ind_up_1] - net_change # perhaps not the underscore2

print('First bit done!')

#Now run a few calculations to determine the fractional change in the adjacent bin 
d1 = np.copy(dens1_plot)
frac_change = np.zeros(len(d1))
for i in range(0, len(d1)):
    frac_change[i] = (net_change[i])/d1[i]
frac_remain = 1.0 - frac_change 
frac_check = (1.0 + rp1**2)/(1.0 + rp2**2)
#frac_remain = frac_remain[:-1]

#WOuld be nice to know the range in the fractional change across the cloud

f_c_range = np.amax(frac_change) - np.amin(frac_change)
#NOw multiply previous dens profile to reclaim new one, as an additional test

brand_new_dens = dens1_plot*frac_remain

print('Second bit done')


#Now determine the mass in the cloud for each separate time, determine what gets
#fed into the cloud from the pellet and where it all goes.

#Full range of cloud first
k = 9
rm1 = np.linspace(rp1,rc1 ,2**k + 1)
dm1 = eps*((1.0 + rp1**2)/(1.0 + rm1**2))
mc1_full = 4.0*np.pi*spint.simps(dm1*rm1**2, rm1)
mc1_full = 4.0*np.pi*spint.romb(dm1*rm1**2, dx = rm1[1] - rm1[0])

mp_1 = (4.0/3.0)*np.pi*(rp1**3)
mp_2 =  (4.0/3.0)*np.pi*(rp2**3)

rm2 = np.linspace(rp2, rc2, 2**k +1)
dm2 = eps*((1.0 + rp2**2)/(1.0 + rm2**2))
mc2_full = 4.0*np.pi*spint.simps(dm2*rm2**2, rm2)
mc2_full = 4.0*np.pi*spint.romb(dm2*rm2**2,dx = rm2[1] - rm2[0])

#The following tests prove true as of 20/09/2019
mass_check1 = mc1_full + mp_1
mass_check2 = mc2_full + mp_2

mass_diff_cloud = mc2_full - mc1_full 
mass_diff_pellet = mp_1 - mp_2

#Now the more robust check on things

dens1 = eps*((1.0 + rp1**2)/(1.0 + r**2))
dens2 = eps*((1.0 + rp2**2)/(1.0 + r**2))

#Now make the test for conservation of mass in each section

mc2_in = 4.0*np.pi*spint.simps(dens2[ind_low_2:ind_low_1]*r[ind_low_2:ind_low_1]**2, r[ind_low_2:ind_low_1])
mc2_out = 4.0*np.pi*spint.simps(dens2[ind_up_1:ind_up_2]*r[ind_up_1:ind_up_2]**2, r[ind_up_1:ind_up_2])
mc2_mid = 4.0*np.pi*spint.simps(dens2[ind_low_1:ind_up_1]*r[ind_low_1:ind_up_1]**2, r[ind_low_1:ind_up_1])
mc2 = mc2_in + mc2_out + mc2_mid 
big_mofo_check = mc2 - mc2_full


dens1 = dens1[ind_low_1:ind_up_1]
dens2 = dens2[ind_low_2:ind_up_2]

#NOw do the mom-conservation across the cloud 

const = dens1[-1]*rcd[t1]
x = len(dens1)
flow = np.zeros(x)
for i in range(0, x):
    flow[i] = const/dens1[i]
#normalised flow speeds calculated 

"Following code is a test on moving a distribution of prticles"

#Now calculate how far the FRACTION of a stopped distribution of electrons moves in delta_t

test = np.zeros(x)

#####################################################################################################
#Gaussian
cent = int(0.5*(ind_up_1 - ind_low_1))
rcent = r[cent]
sig = 2.5
test[:] = (1.0/(np.sqrt(2.0*np.pi*sig**2.0)))*np.exp(-((r[ind_low_1:ind_up_1]- rcent)**2)/(2.0*sig**2))
"""
"""
######################################################################################################
#Polynomial
#test[:] = 2*r_plot1[:]**2 - r_plot1[:] + 3

######################################################################################################
#Linear 
#test[:] = 0.5*r_plot1[:]
#test density is created - now move it by flow*delta_t

new_pos = flow[:]*delta_t/(len(t)) + r[ind_low_1:ind_up_1] 
new_pos_dens = frac_change*test[:] 
#that's the shifted particles calculated 

#Now to calculate the new density profile, the addition of the remaining particles to shifted ones
#on a common grid.

#First thing to do is to put new points on the common grid.
total_new_dens = np.zeros(ind_up_2 - ind_low_2)
test_old = np.zeros(len(r_plot2))
test_old[ind_low_1 - ind_low_2:ind_up_1 - ind_up_2] = test[:]*frac_remain 
total_new_dens[:] += test_old[:]
g = spit.interp1d(new_pos, new_pos_dens, kind = 'cubic')

thing = (np.abs(new_pos[0] - r[:])).argmin() + 1
npd = g(r[thing :ind_up_2-1])

total_new_dens[-len(npd):] += npd[:]

#Test to confirm that nothing from the test profile has been lost

check1 = spint.simps(test,r_plot1, dx = dr)
check2 = spint.simps(total_new_dens, r[ind_low_2:ind_up_2], dx = dr)

"""
for i in range(0,x):
    temp= (np.abs(new_pos[i] - r[:])).argmin()
    new_pos[i] = r[temp]
listy = []
#Also need to interpolate shifted densities from existing grid onto all points of common grid to avoid sawtoothing
g = spit.interp1d(new_pos,new_pos_dens,kind = 'cubic')
#real_new_pos_dens = 

total_new_dens[ind_low_1 - ind_low_2:ind_up_1 - ind_up_2] = dens_old[ind_low_1 - ind_low_2:ind_up_1 - ind_up_2]*frac_remain[:] #+ test[i]*frac_change[i]
for i in range(0,x):
    temp2 = (np.abs(new_pos[i] - r[ind_low_2:ind_up_2])).argmin() 
    total_new_dens[temp2] += new_pos_dens[i]
"""

##########################################################################################

print('Mass check done!')
#Plotting stuff!
fig, ax = plt.subplots() 
ax2 = ax.twinx()

#Vertical line markers to track pellet and cloud positions
ax.axvline(r[ind_low_1], linestyle = 'dotted')
plt.text(0.05,0.14 ,r'$\tilde{r}_p(\tilde{t}_2)$', fontsize = 6)

ax.axvline(r[ind_low_2], linestyle = 'dotted')
plt.text(1.1,0.16, r'$\tilde{r}_p(\tilde{t}_1)$', fontsize = 6)

ax.axvline(r[ind_up_2], linestyle = 'dotted')
plt.text(14.3, 0.03, r'$\tilde{r}_c(\tilde{t}_2)$', fontsize = 6)
#ax.plot(r_plot1, net_change, label = 'change', color = 'firebrick')
#ax.scatter(r_plot1, brand_new_dens, label = 'br_dens', color = 'magenta', marker = 'x')
#ax.scatter(r_plot1,dens1_plot, label = 'first', marker = 'x')
#ax.scatter(r_plot2, dens2_plot, label = 'second', marker = '|')
#ax.scatter(r_plot1, new_dens, label = 'new', marker = '_')
l3 = ax2.plot(r_plot1, test, color = 'black', linestyle = '--', label = r'$n_e(\tilde{t}_1)$')
l1 = ax.plot(r_plot1,dens1_plot, label = r'$\tilde{\rho}_c(\tilde{t}_1)$')
l2 = ax.plot(r_plot2, dens2_plot, label = r'$\tilde{\rho}_c(\tilde{t}_2)$')
l4 = ax2.plot(new_pos, new_pos_dens, label = r'f$n_e (\tilde{t}_1 \rightarrow \tilde{t}_2)$', linestyle = '--', color = 'gold')
l5 = ax2.plot(r_plot2, total_new_dens, label = r'$n_e(\tilde{t}_2)$',linestyle =  '--', color = 'blue')

lines = l1 + l2 + l3 + l4 + l5
labs = [l.get_label()  for l in lines]
ax.set_yscale('log')
#ax2.set_yscale('log')
#ax.set_xscale('log')
ax.set_ylabel(r'$\tilde{\rho}_c$', rotation = 0, fontsize = 14)
ax.set_xlabel(r'$\tilde{r}$', fontsize = 14)
ax2.set_ylabel(r'$\tilde{n}_e$', fontsize = 14, rotation = 0)

ax.yaxis.set_label_coords(-0.06, 0.41)
ax2.yaxis.set_label_coords(1.05, 0.57)
#ax.set_title('Cloud density changes between times '+  r'$\tilde{t} = $'+ str(round(t1/len(t),3))+ ' and ' + r'$\tilde{t} = $' + str(round(t2/len(t),3)))
plt.legend(lines, labs)
plt.grid(which = 'both', axis = 'both')
plt.savefig(mydir+'cloud_flow_test.png', format = 'png', dpi = 1600)
plt.show()