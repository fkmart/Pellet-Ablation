import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rp, rc, eps ,t, dt
import math as mt

t1 = 40
t2 = 60

shift = t2-t1

rp1 = rp[t1] 
rp2 = rp[t2]
rc1 = rc[t1]
rc2 = rc[t2]

r1 = np.linspace(rp1,rc1, num = 500, endpoint = 'true')
r2 = np.linspace(rp2, rc2, num = 500, endpoint = 'true')
dens1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)
dens2 = eps*(1.0 + rp2**2)/(1.0 + r2**2)
dens2_1 = eps*(1.0 + rp2**2)/(1.0 + r1**2)

rpd1 = -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1.0)))))*(mt.tanh(1.0 - t[t1]))/(np.sqrt(np.log(mt.cosh(1.0 - t[t1]))))

dpdt = 2.0*eps*rp1*rpd1/(1.0 + r1**2)

d_dens = dpdt*shift*dt

frac_av = (dens1 + d_dens)/dens1 
print(frac_av[0])
dens_out = np.copy(dens1)

frac_eq = dens2_1/dens1 
print('Equilibrium comparative fractions = ' + str(frac_eq[0]))
frac1 = 1.0
for i in range(0,shift):
    rpd =  -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1.0)))))*(mt.tanh(1.0 -
     t[t1 + i]))/(np.sqrt(np.log(mt.cosh(1.0 - t[t1 + i]))))
    r = np.linspace(rp[t1+i], rc[t1+i], num = 500, endpoint = 'true')
    dens = eps*(1.0 + rp[t1 + i]**2)/(1.0 + r1**2)
    dpdt = 2.0*eps*rp[t1 + i]*rpd/(1.0 + r1**2)
    d_dens = dpdt*dt
    frac = (dens + d_dens)/dens
    dens_out = frac*dens_out
    
    frac1 *= frac

print(frac1[0])
fig, ax = plt.subplots()
ax.plot(r1,dens1, label = r'$\tilde{\rho}_1$')
ax.plot(r1,dens2_1, label = r'$\tilde{\rho}_2$')
ax.plot(r1,dens_out,linestyle = '--', label = r'$\tilde{\rho}_{1 \rightarrow 2}$')
ax.plot(r1, frac_av*dens1, linestyle = '--', label = r'$\tilde{\rho}$')
plt.legend()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.set_ylabel(r'$\tilde{\rho}$', fontsize = 12, rotation = 0)
#ax.set_xlim(left = 0.9, right = 1.5)
#plt.ylim(bottom = 0.006, top = 0.01)
#plt.yscale('log')
#plt.savefig('revised_transport_zoom.png', format = 'png', dpi = 1200)
plt.show()