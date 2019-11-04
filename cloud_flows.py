import numpy as np
from stop_calc_rp_rc import *
import matplotlib.pyplot as plt
import scipy.interpolate as spint

rp = calc_rp(t)
print('rp done')
rc = calc_rc(t)
print('rc done')
rcd, r_i = rdot(t,rc,rp)
print('rcd done')

"""Now time to compare the differences between the t1,t2 and t1+dt density profiles"""

t1 = 100
t2 = 105
dt = t[t2] - t[t1]

def dens(t,r_i):
    den = eps*((1.0 + r_i[t,0]**2)/ (1.0 + r_i[t,:]**2))
    return den
d1 = dens(t1, r_i)
d2 = dens(t2,r_i)

r_change = rcd[t1,:]*dt + r_i[t1,:]

fig, ax  = plt.subplots()

print('Distance covered by cloud at pellet surface in dt is ' + str(rcd[t1,0]*dt))
print('Radius of pellet lost due to ablation ' + str(r_i[t1,0] - r_i[t2,0]))

"""
ax.plot(r_i[t2,:], d2, label = r'$t_2$')
ax.plot(r_i[t1,:], d1, label = r'$t_1$')
ax.plot(r_change,d1, label = r'$t_1 + \Delta t$')
ax.set_yscale('log')
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel(r'$\tilde{\rho}$', rotation = 0)
ax.set_title(r'$t_1 = $' + str(t1) + r'$\mu $s ,'+ r'$t_2 = $' + str(t2) + r'$\mu$s,' + r'$\Delta t = $' + str(t2-t1) + r'$\mu$s' )
plt.legend()
plt.show()"""


"""TO determine the difference between one prodile and the advanced profile for common point"""
for i in range(0, 500):
    if (r_i[t1, -1] < r_i[t2,i]):
        r_up = r_i[t2, i-1]
        i_up = i-1
        break
    else:
        pass

r_low = r_i[t2,1]
i_low = 1

f = spint.interp1d(r_i[t1,:], d1[:], kind = 'cubic')

r_int = np.linspace(r_low, r_up, 500 - (500 - i_up) - i_low)
dens_int = f(r_int)
"""
dens_diff = dens_int - d2[i_low : i_up]
ax.plot(r_int, dens_diff)
plt.show()
"""

tp = 900
ax.plot(r_i[tp,:], rcd[tp,:] )
plt.show()
