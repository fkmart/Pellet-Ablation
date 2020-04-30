import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rp, rc, eps ,t, dt
import math as mt
import romberg as ro
import scipy.interpolate as spint 

t1 = 40
t2 = 45

shift = t2-t1

rp1 = rp[t1] 
rp2 = rp[t2]
rc1 = rc[t1]
rc2 = rc[t2]

r1 = np.linspace(rp1,rc1, num = 513, endpoint = 'true')
r2 = np.linspace(rp2, rc2, num = 513, endpoint = 'true')
dens1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)
dens2 = eps*(1.0 + rp2**2)/(1.0 + r2**2)

"Fractional changes at each end"
rc_frac = rc2/rc1 
rp_frac = rp2/rp1 

"Mid point of r1 array"
mid_ind = int(0.5*(len(r1)))

r1_mid = r1[mid_ind]
r_new = np.copy(r1)
r_new[mid_ind+1:] *= rc_frac
r_new[:mid_ind] *= rp_frac

"Interpolation r array"
r_int = np.linspace(rp2, rc2, num = 1025, endpoint= 'true')

"#Various test mass profiles"
#test_mass = np.zeros(len(r1))
test_mass = (1e-3)*np.exp(-((r1 - r1_mid)/2)**2)
l1 = 100

#test_mass[mid_ind - int(0.5*l1): mid_ind + int(0.5*l1)] = np.ones(l1)*1e-4

"Initial integration and renormalisation for first attempt"
test_integ = ro.romberg_samp(test_mass,r1)
trans_integ = ro.romberg_samp(test_mass,r_new)
renorm = test_integ/trans_integ

h = spint.interp1d(r_new,renorm*test_mass,kind = 'cubic', fill_value = 'true')
test_trans = h(r_int)

"Ratio of elec to neutral density at common points prior to transport"
frac_pre = test_mass/dens1

# now need to interpoate these fractions to the tranported co-ordinates and 
#extrapolate to new cloud and pellet regions
g = spint.interp1d(r1,frac_pre, kind = 'cubic', fill_value = 'extrapolate')
frac_trans = g(r_int)
dens_trans = eps*(1.0 + rp2**2)/(1.0 + r_int**2)
test_trans2 = frac_trans*dens_trans

#interpolate the transformed stuff to maintain good resolution 

#test_pchip = spint.pchip_interpolate(r_new,test_trans, r_int)

#conservation check 
integral1 = ro.romberg_samp(test_mass,r1)
integral2 = ro.romberg_samp(test_trans2,r_int)
integral3 = ro.romberg_samp(test_trans,r_int)
print(integral1)
print(integral2)
print(integral3)
fig,ax = plt.subplots(figsize = (13,8))
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel(r'$\tilde{\rho}_c$', rotation = 0)
ax2 = ax.twinx() 
line = [] 
ax2.set_ylabel(r'$\tilde{n}_s$', rotation = 0)
line += ax.plot(r2,dens2, label = r'$\tilde{\rho}_c(\tilde{t}_2)$')
line += ax.plot(r1,dens1, label = r'$\tilde{\rho}_c (\tilde{t}_1)$')
line += ax2.plot(r1,test_mass, color = 'navy', label = 'Test Seeded Density - ' + r'$\tilde{\rho}_{s,0}$')
line += ax2.plot(r_int, test_trans2, color = 'limegreen', label = 'Fraction conserved transport - ' + r'$\tilde{\rho}_{s,f}$')
line += ax2.plot(r_int, test_trans, color = 'firebrick', label = 'Charge conserved transport - ' + r'$\tilde{\rho}_{s,c}$')
labels = [l.get_label() for l in line]
#plt.yscale('log')
plt.legend(line,labels)
plt.text(1.0,0.002, r'$\leftarrow$' + 'To pellet')
plt.text(5,0.00075, r'$\int_{\tilde{r}_p}^{\tilde{r}_c} \tilde{\rho}_{s,0} = %.8f$' % integral1)
plt.text(5,0.000675, r'$\int_{\tilde{r}_p}^{\tilde{r}_c} \tilde{\rho}_{s,f} = %.8f$' % integral2)
plt.text(5,0.00060, r'$\int_{\tilde{r}_p}^{\tilde{r}_c} \tilde{\rho}_{s,c} = %.8f$' % integral3)
plt.savefig('stretchy_transport.png', format = 'png', dpi = 1600)
plt.show()