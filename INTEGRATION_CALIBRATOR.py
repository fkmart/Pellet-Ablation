import numpy as np
import scipy.integrate as spint

#import stop_calc_rp_rc
from gen_var import * 

tchoice = 160
r0n = 1.0
M0 = (4.0*np.pi *r0n **3)/3.0


MP = (4.0/3.0)*np.pi*rp[tchoice]**3 


l = 500
r = np.linspace(rp[tchoice], rc[tchoice], l)
rho_c = eps*((1.0 + rp[tchoice]**2)/(1.0 + (r[:])**2))
MC = (4.0)*np.pi*spint.simps(rho_c[:]*r[:]**2, r)
mass = lambda x: eps*((1.0 + rp[tchoice]**2)/(1.0 + x**2))*x**2

print('The original total mass is : ' + str(M0))
print('The current mass in the pellet is : ' +str(MP))
print('Mass lost from pellet is ' + str(M0 - MP))

print('-------------------------------------------------')
print('Mass in cloud is ' + str(MC))

#Try look at end times

r_end = np.linspace(rp[-1], rc[-1], l)
rho_end = eps*((1.0 + rp[-1]**2)/(1.0 + (r_end)**2))

MC_end = 4.0*np.pi*spint.simps(rho_end[:]*r_end[:]**2, r_end)
print('The mass in the cloud at the end is : ' + str(MC_end))
