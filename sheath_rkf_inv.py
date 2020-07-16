import numpy as np 
import matplotlib.pyplot as plt 
import rkf45_sheath_inv as rkf 
import scipy.interpolate as spint

"""Real units"""
e = 1.6*10**(-19) # in C
kb = 1.38*10**(-23) # in J/K
mi = 1.67*10**(-27) #in kg
n0 = 10**19 #in m^-3
eps0 = 8.854*10**(-12) # in F/m
me = 9.11*10**(-31) # in kg
phi_0 = 1000.0 # in V


T = phi_0*e/kb # in K
u0 = np.sqrt(kb*T/mi) # in m/s
ld = (eps0*kb*T/(e*e*n0))**(0.5) # in m
epsp = mi*u0*u0/(2.0*e) #in V


"""Normalised Units"""

T_n = T*kb/e #in V
T_n /= phi_0 #non-dim
epsp_n = 20.0*epsp/phi_0 #non-dim
e_n = e/e #non-dim
n_n = n0/n0 # non-dim

#FInal constant
K = (2.0*n0*e*ld*ld/(phi_0*eps0))**(-0.5)

"""Setting up the function to be integrated"""

def sheath_final(phi_n,x):
    term1 = T_n*np.exp(phi_n/T_n) - T_n
    term2 = 2.0*epsp_n*(1.0 - phi_n/epsp_n)**(0.5) - 2.0*epsp_n
    diff1 = term1 + term2 
    diff = K*(diff1)**(-0.5) # produces dx/dphi
    return diff

"""Setting up the integrand matrix"""
phi_s = np.log(np.sqrt(2.0*mi/(2.0*np.pi*me)))
print(phi_s)
#l = 3000
phi_1 = 1e-3
#phi = -np.linspace(phi_s, phi_1, l)
#phi = np.logspace( np.log10(phi_s),np.log10(phi_1), l)
x0 = 0.0
h = -0.01

#epsp_n *=phi_s*1.01

phi, x, dr_arr, err_arr = [],[],[],[]
phi = [0.0]
x = [x0]
deriv_kick = 1e-1
phi_start = deriv_kick*h
phi = np.append(phi,phi_start)
x = np.append(x,h)
"""Normalised integration"""
while np.abs(phi[-1]) < np.abs(phi_s):
    phi_out, x_out, dr_out, err_out = rkf.rkf(phi[-1],x[-1],h,sheath_final,-phi_1,-phi_s)
    phi = np.append(phi, phi_out)
    x = np.append(x, x_out)
    dr_arr = np.append(dr_arr, dr_out)
    err_arr = np.append(err_arr, err_out)

#x *= ld*K
#phi *= phi_0

#deriv = (sheath_final(phi[-1],x[-1]))**-1
#phi_new = np.linspace(phi[-1], 0.5,num = 500)
#x_new = spint.pchip_interpolate(phi,x,phi_new)

fig, ax = plt.subplots()
ax = plt.plot(x,np.abs(phi),color = 'navy')
#plt.yscale('log')
plt.show()