import numpy as np 
from gen_var import rp,rc, eps, t
import matplotlib.pyplot as plt 
import math as mt 
from stop_calc_rp_rc import u_solution


"""Original method was to move particles forward by 1 step and an arbitrary time
 - this was calculated by rp/rpdot for the selected time. """

"""Revised method - determine u and move forward by the distance corresponding to dt*u"""

t1 = 250 
t2 = 270

rp1, rc1 = rp[t1], rc[t1]
rp2, rc2 = rp[t2], rc[t2] 

shift = t2 -t1 
dt = t[1] - t[0]

tshift = dt*shift

l = 400
r1 = np.linspace(1.0, rc1, num = l)
dr = r1[1] - r1[0]
r2 = np.append(r1,rc1 + dr)
rho1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)

rho2_theory = eps*(1.0 + rp2**2)/(1.0 + r1**2)

rpdot1 = -(0.5*mt.tanh(1.0 - t[t1]))/(np.sqrt(np.log(mt.cosh(1)))*np.sqrt(np.log(mt.cosh(1 - t[t1]))))

a_time = rp1/rpdot1
print(tshift)
drhodt = -2*eps*rp1*rpdot1/(1.0 + r1**2)

rho2_emp = -drhodt*np.abs(tshift) + rho1

#fractional change in neutral density 

drho1 = drhodt*tshift 
frac_lost = drho1/rho1 
frac_rem = 1 -frac_lost

#define electron distribution 

#gaussian shape 
A = 0.005 
sig = 0.625 
r_cent = 3.0 

e_rho = A*np.exp(-(r1 - r_cent)**2/sig**2)

#Need to now add the space moved by the cloud to co-ordinates
u = u_solution(t1,t2,l)
x_move= dt*shift*u

r_trans = r1 + x_move

e_rho_trans = np.zeros(len(r2))
e_rho_trans[0] = frac_rem[0]*e_rho[0]

"""New method""" 

dens_rem = np.zeros(l )
dens_lost = np.zeros(l) 

dens_rem[:] = frac_rem[:] * e_rho 
dens_lost[:] = frac_lost[:] * e_rho 

e_trans = A*np.exp(-(r_trans - r_cent)**2/sig**2)
e_trans *= frac_rem[0]
e_trans += dens_lost

"""" - The odl method
for i in range(1, len(rho1)):
    e_rho_trans[i] = frac_rem[i]*e_rho[i] + frac_lost[i]*e_rho[i-1]
e_rho_trans[-1] = -frac_lost[i]*e_rho[-1]
r_trans = np.append([r1[0]], r_trans)"""

print('Now test to check that summed along r is conserved')
print('Initial sum of electrons ' + str(np.sum(e_rho)))
print('Post-transport sum of electrons ' + str(np.sum(e_rho_trans)))
print('t1 density at rp is ' + str(rho1[0]))
print('Change in density at rp ' + str(drho1[0]))
print('Fraction lost which is a constant ' + str(frac_lost[0]))
print('Fraction remaining which is a constant ' + str(frac_rem[0]))
fig,ax = plt.subplots()
ax.plot(r1, rho1, label = r'$\tilde{\rho}(\tilde{t}_1)$')
ax.plot(r1, rho2_emp, label = r'$\tilde{\rho}(\tilde{t}_2)$' + ' - empirical')
ax.plot(r1, e_rho, label = r'$\tilde{n}_e(\tilde{t}_1)$')
ax.plot(r_trans, e_trans, label = r"$\tilde{n}_e(\tilde{t}_2)$" + ' - transported')
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel('Normalised Particle Density')
plt.legend()
plt.text(20,0.004, r'$\rightarrow$' + 'to plasma')
#plt.savefig('excel_transport.png', format = 'png', dpi = 1400)
plt.show()