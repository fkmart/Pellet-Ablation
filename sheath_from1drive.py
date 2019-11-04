import numpy as np
from scipy import integrate as spit
import matplotlib.pyplot as plt 

Te = 10.0**(-1) #Real - actually kb*T/e
epsp = 1.00 #Non-dim
n = 10**20 #Real
esp = 8.854*10**(-12) #Real in F/m
me = 9.11*10**(-31) #Real in kg
mp = 1.67*10**(-27) #Real in kg
phiw = np.log((mp/(2.0*np.pi*me))**(0.5)) #non-dim
dp = 0.001
 
phis = -np.logspace(0.0001, phiw , 500) #Non-dim

kick = 1e-6
points = 5000
phis = -np.logspace(np.log10(kick), np.log10(phiw),points ) #non-dim
#phis = np.flip(phis, axis = 0)
e = 1.6*10**(-19) #Real 

#epsp = -phis[-1]*1.05
def dxdp(x,phi):
    #epsp = 1.0 - phi
    dxdphi = (np.exp(phi) - 1.0 + 2.0*epsp*(1.0 - phi/epsp)**(0.5) - 2.0*epsp)**(-0.5)
    return dxdphi

cofactor = (2*e*n/(esp*Te))**(-0.5)

k1,k2,k3,k4 = 0.0, 0.0, 0.0, 0.0
xn0 = 0.0
x_out = []
x_out.append(xn0)
for i in range(0, len(phis) - 1):
    xn0 = x_out[i]    
    dp = phis[i+1] - phis[i]
    k1 = dp*dxdp(xn0, phis[i])
    k2 = dp*dxdp(xn0 + 0.5*dp, phis[i] + 0.5*k1)
    k3 = dp*dxdp(xn0 + 0.5*dp, phis[i] + 0.5*k2)
    k4 = dp*dxdp(xn0 + dp, phis[i] + k3)
    xn1 = xn0 + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
    x_out.append(xn1)
x_out = np.asarray(x_out)*cofactor 
ld = np.sqrt(Te*esp/(e*n))

x = x_out[1:]/ld
x = x_out/ld
plt.xlabel(r'$\frac{x}{\lambda_D}$', fontsize = 14)
plt.ylabel(r'$\frac{\phi_s}{T_e}$', fontsize = 14, rotation = 0)
#plt.xscale('log')
#plt.yscale('log')

ax = plt.gca()
ax.xaxis.set_label_coords(0.55, -0.055)
ax.yaxis.set_label_coords(-0.08, 0.50)
plt.plot(x,phis)
plt.show()
#plt.savefig('Sheath_pot.png', format = 'png', dpi = 1600)

print('Sheath size is', x[-1] - x[0])

s = spit.simps(dxdp(0.0,phis), phis)
print(s)