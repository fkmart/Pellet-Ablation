import numpy as np
import iterative_sol 
from gen_var import t_static, dt, t, F0, e, epsilon0, dr, pel_pot, pel_dens_numb, r0
from electron import RME, M_fac
import stop_calc_rp_rc 
import discret_mat
import interpolate_limits
import scipy.integrate as spint
import matplotlib.pyplot as plt
import os

rp = stop_calc_rp_rc.calc_rp(t)
rp_static = rp[t_static]
rc = stop_calc_rp_rc.calc_rc(t)
rc_static = rc[t_static]

number = F0*np.pi*dt*(rp_static**2)*10**(-9)
#assume that all occurs inside a cylinder of radius rp_static
real_e_dens = number/((np.pi*(rp_static*10**(-3))**2)*(rc_static-rp_static)*10**(-3)) # as a number density

"can now make a calculation of the potential field in this domain"
mydir = './static_outputs_phi' 
mydir = os.path.join(os.getcwd(), 'static_outputs_phi') + os.sep
p = 0
d_file = np.loadtxt(mydir+'/real_density_t'+str(t_static) +'pot'+str(pel_pot[p])+'.txt')

dens, cloud_r = interpolate_limits.interpol(d_file)

"""Some tests to be performed first
1) Sum the inteprolated density via simpsons acorss cloud - check
this equals the value found by the same process on the original
2)Calculate the actual number density from number of electrons
in cylinder and integrate via simpsons to confirm fractino of total 
number is equal to previous calculation"""

check1 = spint.simps(d_file[:,0],d_file[:,1])
check2 = spint.simps(dens,cloud_r)

#Now try to determine the actual number density
actual_dens = dens*number 

#do a conservation test here
check3 = spint.simps(actual_dens,cloud_r)
check4 = check2*number

actual_dens /= (dr*np.pi*rp_static**2)
actual_dens /= 10**(-9) 

phi = np.zeros(len(cloud_r))
A = discret_mat.discret(cloud_r)
K = 1.0/(pel_dens_numb*RME*M_fac)
inp = K*actual_dens/epsilon0
#now to normalise RHS of equation so that output is correct.
phi = iterative_sol.SOR(A,phi,inp ,cloud_r)

#Now we can test to see if the solver is right or not
check_again = np.dot(A,phi)
check_again /= dr**2

fig, ax  = plt.subplots()
ax.plot(cloud_r, phi, label = 'potential')
ax2 = ax.twinx() 
ax2.plot(cloud_r, inp)
plt.show()
print('Did we do it?')