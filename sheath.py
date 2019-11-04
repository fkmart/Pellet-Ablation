import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt
import iterative_sol as it
import discret_mat as dm 
import RK4

"""Real units"""
e = 1.6*10**(-19) # electronic charge in C
kb = 1.38*10**(-23) # boltzmann constantin J/K
mi = 1.67*10**(-27) #proton mass in kg
n0 = 10**14 #plasma density in m^-3
eps0 = 8.854*10**(-12) # vacuum permittivity in F/m
me = 9.11*10**(-31) #electron mass  in kg
phi_0 = 10**3 # Plasma temperature in V


T = phi_0*e/kb # Temperature in K
u0 = np.sqrt(2*kb*T/mi) # initial particle velocity in m/s
ld = (eps0*kb*T/(e*e*n0))**(0.5) # Debye length in m
epsp = mi*u0*u0/(2.0*e) # initial particle energy in V


"""Normalised Units"""

T_n = T*kb/e # Temperature in V (again)
T_n /= phi_0 #non-dim
epsp_n = epsp/phi_0 #non-dim
e_n = e/e #non-dim
n_n = n0/n0 # non-dim


#FInal constant
K = (2.0*n0*e/(phi_0*eps0))**(-0.5)

"""Setting up the function to be integrated"""

def sheath_final(phi_n):
    term1 = T_n*np.exp(phi_n/T_n) - T_n
    term2 = 2.0*epsp_n*(1.0 - phi_n/epsp_n)**(0.5) - 2.0*epsp_n
    diff1 = term1 + term2 
    diff = (diff1)**(-0.5)
    return diff

def sheath_real(phi):
    term1 = phi_0*np.exp(phi/phi_0) - phi_0
    term2 = 2.0*epsp*(1.0 - phi/epsp)**(0.5) - 2.0*epsp
    diff1 = term1 + term2
    diff2 = 2.0*diff1*e*n0/eps0
    diff3 = diff2**(-0.5)
    return diff3


"""Setting up the integrand matrix"""
phi_s = np.log(np.sqrt(mi/(2.0*np.pi*me)))
print(phi_s)
l = 5000
phi_1 = 1e-6
phi = -np.linspace(phi_1, phi_s, l)
phi = -np.logspace(np.log10(phi_1), np.log10(phi_s), l)
#phi = np.flip(phi, axis = 0)
x0 = 0.001


epsp_n *=-phi_s*1.05
epsp_n = 1.00
epsp *= -phi_s*1.05

"""Normalised integration"""
x = RK4.RK4_y(sheath_final, x0, phi) 
#x *= K
#phi *= phi_0

fig, ax = plt.subplots()
ax.plot(x, phi,label = 'norm')
ax.set_xlabel(r'$\frac{x}{\lambda_D}$', rotation = 0, fontsize = 14)
ax.set_ylabel(r'$\frac{\phi}{T_e}$', rotation = 0, fontsize  =14)
ax.xaxis.set_label_coords(0.53,-0.05)
ax.yaxis.set_label_coords(-0.08, 0.48)

#ax.plot(x_real*(2.0*e*n0/eps0)**(-0.5) /ld, phi_real, label = 'real')
#ax.set_yscale('log')
#ax.set_xscale('log')
print('Sheath size is ' + str(x[0] - x[-1]))
plt.grid()
plt.legend()
plt.show()