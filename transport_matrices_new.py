import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rc, rp, t,dt  , eps 
from stop_calc_rp_rc import rdot, rcdot 
from transport_functions import transport
import scipy.interpolate as spint
import find_nearest as fn 

t1 = 40 
t2 = 45 
shift = t2-t1 

rp1 = rp[t1] 
rp2 = rp[t2]
rc1 = rc[t1] 
rc2 = rc[t2] 

l = 1000

r1 = np.linspace(rp1,rc1, num = l, endpoint = 'true')
r2 = np.linspace(rp2, rc2, num = l, endpoint = 'true')

grid = np.linspace(rp2,rc2,num = 10000, endpoint = 'true')

dens1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)
dens2 = eps*(1.0 + rp2**2)/(1.0 + r2**2)

"Now calculate the fluid velocities"

rdots = rdot(t,rc,rp)
rcd = rdots[t1]

const = dens1[-1]*rcd

u1 = const/dens1[:]
trans_r = r1 + dt*shift*u1

"Define matrices for the method"

A = np.zeros(3*l) # multiplicative matrix
A = np.reshape(A,(3,l))
A[2,:]  = 1.0 

dens_comp = np.zeros(3*l) #compsite density matrix
dens_comp = np.reshape(dens_comp, (l,3))

dens1_2 = eps*(1.0 + rp1**2)/(1.0 + r2**2) # evaluate density1 at r2 points and remove outliers
for i in range(0, l):
    if r2[i] < rp1 or r2[i] > rc1:
        dens1_2[i] = 0.0

k = 5.0
dens_comp[:,2] = eps*np.exp(-k*(r2 - rp2)**2)
dens_comp[:,0] = dens1_2[:]

g = spint.interp1d(trans_r,r1,kind = 'cubic')
il = next(p[0] for p in enumerate(r2) if p[1] > trans_r[0])
pre_r1 = g(r2[il:])
pre_r = np.zeros(l)
pre_r[il:] = pre_r1
pre_dens = eps*(1.0 + rp1**2)/(1.0 + pre_r**2)

for i in range(0, len(pre_r)):
    if pre_r[i] < rp1 or pre_r[i] > rc1:
        pre_dens[i] = 0.0

dens_comp[:,1] = pre_dens[:]
f_rem = np.zeros(l)
f_trans = np.zeros(l)

for i in range(l-1,il,-1): 
    f_trans[i] = (dens2[i] - np.dot(A[:,i],dens_comp[i,:]))/pre_dens[i]
    rem = 1.0 - f_trans[i] 
    k = fn.find_nearest(r2,pre_r[i])
    A[0,k] = rem 
    A[1,i] = f_trans[i]

dens_new = np.zeros(l)  
for s in range(0,l):
    dens_new[s] = np.dot(A[:,s], dens_comp[s,:])
fig,ax = plt.subplots() 
ax.plot(r2,dens2)
ax.plot(r2, dens_new)
plt.yscale('log')
plt.show()

fig,ax = plt.subplots() 
ax.plot(pre_r,A[0,:])
plt.show()