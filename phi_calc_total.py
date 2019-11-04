import elec_phi_calc
import numpy as np
import charge_dens
import scipy.integrate as spinteg
import matplotlib.pyplot as plt
from electron import RME, M_fac
from gen_var import dr, t
import os
import time

def phi_tot(t):
    tick1 = time.time()
    q_dens, r, norm_const = charge_dens.e_dens(t)
    print(spinteg.simps(q_dens, r))
    phi_wall = 0.0
    phi_plas = 0.0
    phi_e,phi_mat_e = elec_phi_calc.elec_phi(r, q_dens, phi_wall, phi_plas)
    phi_e *= -1.0
    phi_mat_e *= -1.0
    depth = 5 #grid points into cloud that ions penetrate - arbitrarily chosen
    dens_i = np.zeros(len(r))
    dens_i[-depth:] = 1.0
    tot = spinteg.simps(dens_i[-depth:],r[-depth:])
    dens_i[-depth:] /= tot
    print(spinteg.simps(dens_i[-depth:], r[-depth:]))
    phi_i, phi_mat_i = elec_phi_calc.elec_phi(r, dens_i, phi_wall, phi_plas)
    phi_i*= -1.0
    phi_mat_i *= -1.0
    #phi_e /= r[-1]**2
    #phi_i /= (depth*dr)**2
    phi_tot = phi_i + phi_e
    #np.savetxt(os.path.join('Analysed Outputs', 'elec_phi_t'+str(t)+'.txt'), phi_e)
    #np.savetxt(os.path.join('Analysed Outputs', 'ion_phi_t'+str(t)+'.txt'), phi_i)
    tock1 = time.time()
    print('Time for outer loop is ' + str(np.abs(tick1 - tock1)))
    return phi_tot, r, norm_const, phi_e, phi_i, q_dens, dens_i, phi_mat_e, phi_mat_i
#for i in range(1, len(t)):
x, r, K, pe, pi,q_den, i_den, pe_mat, pi_mat = phi_tot(1)


x = np.append(x,0.0)
pe = np.append(pe, 0.0)
pi = np.append(pi, 0.0)
r = np.append(r, r[-1] + dr)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ln1 = ax.plot(r,x*RME*M_fac*K, '--', label = r'$\phi_{t,1}$')
#ln2 = ax.plot(r, pi*RME*M_fac*K, label = 'ions1')
#ln3 = ax.plot(r, pe*RME*M_fac*K, label = 'electron1')
#ln4 = ax.plot(r[1:], (pe_mat+pi_mat)*RME*M_fac*K, label = 'matrix1')
ax.set_ylabel('$\phi$/V', fontsize = 13, rotation = 0, x= -0.05, y = 0.57)
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, x = 0.55, y = -0.2)
#ax.set_yscale('log')

ax2 = ax.twinx()
ln9 = ax2.plot(r[1:], (i_den + q_den),'-', label = r'$\rho_{t,1}$')
ln10 = ax2.plot(r[1:], i_den, '-' , label = r'$\rho_{i,1}$')
ln11 = ax2.plot(r[1:], q_den, '-', label = r'$\rho_{e,1}$')

#plt.plot(r[1:], pe_mat*RME*M_fac*K, label = 'e matrix')
#plt.plot(r[1:], pi_mat*RME*M_fac*K, label = 'i matrix')
#plt.plot(r[1:], (pi_mat+pe_mat)*RME*M_fac*K, label = 'total matrix')
#plt.yscale('log')

x2, r2, K2, pe2, pi2,q_den2, i_den2, pe_mat2, pi_mat2 = phi_tot(2)
r2 = np.append(r2, r2[-1] + dr)
r2 = r2[1:]
ln5 = ax.plot(r2,x2*RME*M_fac*K2, '--',  label = r'$\phi_{t,2}$')
#ln6 = ax.plot(r2, pi2*RME*M_fac*K2, '--', label = 'ions2')
#ln7 = ax.plot(r2, pe2*RME*M_fac*K2,'--', label = 'electron2')
#ln8 = ax.plot(r2, (pe_mat2+pi_mat2)*RME*M_fac*K2,'--', label = 'matrix2')


ln12 = ax2.plot(r2, (i_den2 + q_den2), '-', label = r'$\rho_{t,2}$')
ln13 = ax2.plot(r2, i_den2, '-', label = r'$\rho_{i,2}$')
ln14 = ax2.plot(r2 , q_den2, '-', label = r'$\rho_{e,2}$')

ax2.set_ylabel(r'$\rho$ (arbitrary units)')

#lns = ln1 + ln2 + ln3 + ln4 + ln5 + ln6 + ln7 + ln8 + ln9+ln10 + ln11 + ln12 +ln13 + ln14
lns = ln9 + ln10 + ln11 + ln12 + ln13 + ln14 + ln1 + ln5
labs = (l.get_label() for l in lns)
ax.legend(lns, labs, loc = 3)
plt.show()
