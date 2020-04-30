import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rc, rp, t,dt  , eps 
from stop_calc_rp_rc import rdot, rcdot 
from transport_functions import transport
import scipy.interpolate as spint

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
dens_comp[:,2] = eps*np.exp(-k*(r2 - rp2))
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
index_low = next(p[0] for p in enumerate(r2) if p[1] > rc1)-1
index_up = l 
diff = index_up - index_low
f_rem = np.zeros(l)
f_trans = np.zeros(l)

trans = np.zeros(diff)
for i in range(index_low,l):
    f_trans[i] = (dens2[i] - np.dot(A[:,i],dens_comp[i,:]))/pre_dens[i]
f_rem[index_low - diff:index_low] = 1.0 - f_trans[-diff:]
#Update A, multiplcation matrix 
A[0,index_low - diff:index_low] = f_rem[index_low-diff:index_low]
A[1,index_low:] = f_trans[index_low:]
pf = [] 
pf = f_trans[index_low:]
pr = pre_r[index_low:]
pd_final = dens2[index_low:]
pd_comp = np.zeros(l)
for i in range(index_low,l):
    pd_comp[i] = np.dot(A[:,i],dens_comp[i,:])
fig, ax = plt.subplots()
#ax.plot(pr,pf)
ax.plot(pr,pd_final)
ax.plot(pr, pd_comp[index_low:])
plt.show()

index_up = np.copy(index_low)
#index_low = next(p[0] for p in enumerate(r2) if p[1] > trans_r[index_up])
index_low = index_up - diff
count = 0
print(il)
print('Constant difference is' + str(diff))
while index_low > il:
    print('Current lower index is' + str(index_low))
    for i in range(index_low,index_up):
        f_trans[i] = (dens2[i] - np.dot(A[:,i],dens_comp[i,:]))/pre_dens[i]
    f_rem[index_low-diff:index_low] = 1.0 - f_trans[index_low:index_up]
    #Update A, multiplcation matrix 
    A[0,index_low-diff:index_low] = f_rem[index_low-diff:index_low]
    A[1,index_low:] = f_trans[index_low:] 
    pf = f_trans[index_low:]
    pr = r2[index_low:]
    pd_final = dens2[index_low:]
    pd_comp = np.zeros(l)
    for i in range(index_low,l):
        pd_comp[i] = np.dot(A[:,i],dens_comp[i,:])
    fig, ax = plt.subplots()
    #ax.plot(pr,pf)
    ax.plot(pr,pd_final)
    ax.plot(pr, pd_comp[index_low:], linestyle = '--')
    plt.show()
    index_up = np.copy(index_low)
    #index_low = next(p[0] for p in enumerate(r2) if p[1] > trans_r[index_up])
    index_low = index_up - diff

"""Once the last jump takes us to one jump away from end of array it is
time to start considering the particular solutions at the ends
"""
diff = index_low - il 
print('Current lower index is' + str(index_low))
for i in range(index_low,index_up):
    f_trans[i] = (dens2[i] - np.dot(A[:,i],dens_comp[i,:]))/pre_dens[i]
#differences no longer match - must think here
f_rem[index_low-diff:index_low] = 1.0 - f_trans[index_low:index_up]
#Update A, multiplcation matrix 
A[0,index_low-diff:index_low] = f_rem[index_low-diff:index_low]
A[1,index_low:] = f_trans[index_low:] 
pf = f_trans[index_low:]
pr = r2[index_low:]
pd_final = dens2[index_low:]
pd_comp = np.zeros(l)
for i in range(index_low,l):
    pd_comp[i] = np.dot(A[:,i],dens_comp[i,:])
fig, ax = plt.subplots()
#ax.plot(pr,pf)
ax.plot(pr,pd_final)
ax.plot(pr, pd_comp[index_low:], linestyle = '--')
plt.show()