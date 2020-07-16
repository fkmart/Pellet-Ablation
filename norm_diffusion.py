import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as scisig 
import romberg as romb

T_n = 11600 
r0c = 0.1 
r0A = 1e7 
pT = 1.0 
M0 = 1.0 
tf = 1e-3
K_prime = pT*r0A**2 *r0c**2 * np.sqrt(M0)/(T_n**(3.0/2.0) * tf)

M1 = 1.0
M2 = 2.0 
p = 1.0 
sigma = 2.968 #in angstrom 
omega = 1.0
T = 11600.0 
K = 0.001858

M1n = M1/M0 
M2n = M2/M0 
pn = p/pT 
sigman = sigma/r0A
Tn = T/T_n
K_n = K/(K_prime)

numer = (K_n*Tn**1.5 *np.sqrt(1.0/M1n + 1.0/M2n))
denom = (pn * sigman**2 * omega)

D_n = (K_n*Tn**1.5 *np.sqrt(1.0/M1n + 1.0/M2n)) /(pn * sigman**2 * omega)

D_real = D_n * (tf/(r0c**2))**(-1)

#Now to carry out the diffusion 

l = 257
x = np.linspace(1.0, 4.0, num = l, endpoint = 'true')
mid = 2.5
#x = np.linspace(-2.0, 2.0,num = l, endpoint = 'true')
n0 = np.zeros(l)
n0[105:145] = 1.00* np.linspace(1.0,2.0, num = 145-105, endpoint = 'true')

D = D_real
n = n0
count = 1
print('Diffusion coefficient in real units is ' + str(D_real) +  'cm^2/s')
fig, ax = plt.subplots()
ax.plot(x,n0, label = r'$\tilde{\rho}_0$')
plt.xlabel(r'$\tilde{x}$', fontsize = 12)
ax.xaxis.set_label_coords(0.45,-0.04)
plt.ylabel(r'$\tilde{n}$', rotation = 0, fontsize = 12)
ax.yaxis.set_label_coords(-0.06, 0.40)
t = 0.000004

ax.set_yticks(np.arange(0.00,2.50, 0.5))
while count < 6:
    B = romb.romberg_samp(n,x)
    print('Constant equals ' + str(B))
    Z = (x - mid)/(np.sqrt(t*D*4.0))
    f = (np.exp(-Z**2))/(np.sqrt(np.pi*D*t))
    n_diff = scisig.convolve(n, f, mode = 'same')
    j = romb.romberg_samp(n_diff,x)
    print('New integral equals ' + str(j))
    print('Fixed integral equals ' + str(romb.romberg_samp(B*n_diff/j,x)))
    n_diff *= (B/j)
    ax.plot(x,n_diff, label = r'$\Delta \tilde{t} = $' + str(count))
    count +=1
    t += t
    #n = n_diff + n0
    
plt.legend() 
plt.savefig('diff_over_t.png', format = 'png', dpi = 1400)
plt.show()