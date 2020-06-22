import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rp,rc,eps , t
from TRANSPORT_FUNC import transport 

t1 = 40 
t2 = 50 

shift = t2 - t1 
rp1, rc1 = rp[t1], rc[t1] 

rp2, rc2 = rp[t2], rc[t2] 

l = 513 

r1 = np.linspace(rp1,rc1, num = l, endpoint = 'true')
r2 = np.linspace(rp2, rc2, num = l , endpoint = 'true')
rho1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)
rho2 = eps*(1.0 + rp2**2)/(1.0 + r2**2)

rho1_trans, r = transport(t1,shift,l,rho1)

#proxy electron distribution 
e_dens = eps*np.exp(-((r1 - 2.0))**2)
e_dens_trans, r = transport(t1,shift,l,e_dens)

fig,ax = plt.subplots() 
ax.set_xlim(0.75,8)
plt.axvline(rp1, color = 'black', linestyle = ':')
plt.text(0.95,5e-4, r'$r_p(t_1)$', fontsize = 10)
plt.axvline(rc1, color = 'black', linestyle = ':')
plt.text(4.3,0.008, r'$r_c(t_1)$', fontsize = 10)
plt.axvline(rp2, color = 'black', linestyle = ':')
plt.text(0.76,3e-5,r'$r_p(t_2)$', fontsize = 10)
plt.axvline(rc2, color = 'black', linestyle = ':')
plt.text(6.1,0.007, r'$r_c (t_2)$', fontsize = 10)
plt.text(1.00, 1e-6, r'$\leftarrow$' +  ' to pellet')
ax.plot(r1, rho1, label = r'$\rho_n(t_1)$')
ax.plot(r2, rho2, label = r'$\rho_n(t_2)$')
ax.plot(r, rho1_trans, label = r'$\rho_n(t_1 \rightarrow t_2)$')
ax.plot(r1,e_dens, label = r'$\rho_e(t_1)$')
ax.plot(r, e_dens_trans, label = r'$\rho_e(t_1 \rightarrow t_2)$')

plt.yscale('log')
plt.xscale('log')
plt.legend()
#plt.savefig('transported_densities_new.png', format = 'png', dpi = 1400)
plt.show()