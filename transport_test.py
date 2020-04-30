import numpy as np 
import matplotlib.pyplot as plt 
import elec_transport_push as mover
from gen_var import t as t_arr
import stop_calc_rp_rc
from gen_var import rp, rc, eps , dt
import grid_pusher as gp
import romberg as ro
import scipy.interpolate as spint 

t = 50
shift = 20

rp1 = rp[t]
rc1 = rc[t]

rp2 = rp[t+shift]
rc2 = rc[t+shift]

rcd, other = stop_calc_rp_rc.rcdot(t_arr,rc,rp)
grid = np.linspace(rp[t], rc[t], num = 1025, endpoint = 'true')
grid2 = np.linspace(rp2, rc2, num = 1025, endpoint = 'true')

n_dens1 = eps*((1.0 + rp1**2)/(1.0 + grid[:]**2))
n_dens2 = eps*((1.0 + rp2**2)/(1.0 + grid2[:]**2))
const = n_dens1[-1]*rcd[t,-1]
rdot = np.zeros(len(grid))
frac_left = (1.0 + rp2**2)/(1.0 + rp1**2)
for i in range(0, len(grid)):
    rdot[i] = const/n_dens1[i]

frac_change = 1 - frac_left 
r_trans = grid + rdot*dt*shift
n_trans = frac_change*n_dens1
n_rem = frac_left*n_dens1

low = next(p[0] for p in enumerate(r_trans) if p[1] > grid[-1])

dens_out = np.zeros(len(grid))
dens_out[:] += n_rem[:]

for i in range(0, low):
    ind = (np.abs(grid[:] - r_trans[i])).argmin()
    dens_out[ind] += n_trans[i]

dens_out = np.append(dens_out, n_trans[low:])
r_out = np.append(grid,r_trans[low:])

"""#Now need to test for mass conservation between the two times - or atleast
#test that the analytical cloud in some way matches the transport version and match them"""

pre_dens = np.copy(n_dens1) # first analytical profile
pre_int = ro.romberg_samp(4*np.pi*grid*grid*pre_dens,grid)

print('Integrated pre-density ' + str(pre_int))

post_dens = np.copy(n_dens2) # second analytical profile
post_int = ro.romberg_samp(4*np.pi*grid2*grid2*post_dens, grid2)

print('Integrated post-density '+ str(post_int))

#need to interpolate the transported particles onto romberg grid 
g = spint.interp1d(r_trans,n_trans,kind = 'cubic')
trans_grid = np.linspace(r_trans[0], r_trans[-1], num = 257, endpoint = 'true')
trans_dens = g(trans_grid)
#Now to integrate
int_move = ro.romberg_samp(4*np.pi*grid*grid*n_trans,grid)
int_trans = ro.romberg_samp(4*np.pi*trans_grid*trans_grid*trans_dens,trans_grid)
int_rem = ro.romberg_samp(4*np.pi*grid*grid*n_rem,grid)

#Now the pellet's input 
pel_input = (4/3)*np.pi*(rp1**3 - rp2**3 )
print('Mass lost by pellet '+str(pel_input))
print('Mass not transported in cloud 1 ' + str(int_rem))
print('Transported bit ' + str(int_move))
print('Mass after transport ' + str(int_trans))
print('Constituent parts of cloud 1 ' + str(int_move + int_rem))
print('Mass difference of clouds ' + str(post_int - pre_int))

#Now calculate bits at the end and front
g1 = spint.interp1d(grid2,n_dens2, kind = 'cubic')
rps = np.linspace(rp2,rp1,num = 129, endpoint = 'true')
rp_dens = g1(rps)

g2 = spint.interp1d(grid2,n_dens2,kind = 'cubic')
rcs = np.linspace(rc1,rc2, num = 129, endpoint = 'linear')
rc_dens = g2(rcs)

int_rps = ro.romberg_samp(4.0*np.pi*rps*rps*rp_dens,rps)
int_rcs = ro.romberg_samp(4.0*np.pi*rcs*rcs*rc_dens,rcs)
print('Constituent parts of second cloud ' + str(int_rps + int_rcs + int_rem))
print('End bits of second cloud '  +str(int_rcs + int_rps))

"""
fig,ax = plt.subplots()
ax.plot(grid,n_dens1, color = 'navy', label = r'$\tilde{\rho}_c(r,t_1)$')
ax.plot(r_trans,n_trans, color = 'skyblue', label = 'trans')
ax.plot(grid, n_rem, color = 'skyblue', linestyle = '--', label = 'rem')
ax.plot(grid2,n_dens2, color = 'gold', label = r'$\tilde{\rho}_c(r,t_2)$')
ax.plot(r_out, dens_out, color = 'orange', label = 'total')
plt.legend()
plt.show()"""

x = ro.romberg_samp(4*np.pi*trans_grid*trans_grid*(int_move/int_trans)*trans_dens,trans_grid)
print(x)

fig, ax = plt.subplots(figsize = (8,6))
ax.plot(grid, n_trans, color = 'blue', label = 'pre-transport')
ax.plot(trans_grid, trans_dens, color = 'red', label = 'post-transport')
ax.plot(trans_grid, (int_move/int_trans)*trans_dens, color = 'green', label = 'post-transport w/ conserved mass')
ax.set_ylabel(r'$\tilde{\rho}_c$', rotation = 0, fontsize = 12)
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
plt.tight_layout()
plt.legend()
plt.savefig('TRANS_DENS.png', format = 'png', dpi = 1400)
plt.show()