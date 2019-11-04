import numpy as np 
import scipy.integrate as spint
from gen_var import rc, rp, dr, r0, epsilon0
from electron import RME, M_fac 
import SOR_lapl_1D
import matplotlib.pyplot as plt
import os 


t = 1
dens_e = np.loadtxt(os.path.join('Analysed Outputs', 'GNU_electron_density_cloud_1keV_fs_t' +str(t) +'.txt'))


r = np.arange(0, rc[t] - rp[t], dr)

dens = np.zeros(len(r))
ind = 50
phi_out = 1000
phi_pel = -2800

phi_0 = RME*M_fac 

phi_out /= phi_0
phi_pel /= phi_0 

dens[-ind:] = 0.1 

nor = spint.simps(dens[-ind:], r[-ind:])
dens /= nor

norm = r0*(rc[t] - rp[t])**2/(phi_0*epsilon0)
c = spint.simps(dens[-ind:], r[-ind:])
print(c)
rel_tol = 1e-6
phi_i = np.zeros(len(r)-2)
phi_e = np.zeros(len(r) -2)
phi_i = SOR_lapl_1D.SOR(phi_i, dens, rel_tol, phi_pel, phi_out)
phi_e = SOR_lapl_1D.SOR(phi_e, dens_e, rel_tol, phi_pel, phi_out)

plt.plot(r,phi_i*phi_0)
plt.plot(r, phi_e*phi_0)
#plt.plot(r,dens)
plt.show()