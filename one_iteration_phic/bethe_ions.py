from gen_var import *
from proton import RME, M_fac
import matplotlib.pyplot as plt
from electron import RME as RME_e
den = 1.0
nondim = (r0cgs/RME)*solid_dens #non-dimensionalisation factor, in cgs to work with formula
           
E = np.arange(50,20000,20)

T = E/(RME*M_fac)

c1= 0.153536 #factor that must be an amalgamation of many constants
            
#B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
B = -2.0*np.log(I/(RME*M_fac)) #another factor from formula
numer = 2.0*T*RME_e/RME
denom = I/(RME*M_fac)
logarg = numer/denom
B = 2.0*np.log(logarg)
#b = np.log(T*RME*RME_e*M_fac/I)
dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form 
dxde_real = dxde_gas*RME/r0cgs

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(E,-dxde_real)
ax.set_xlabel(r'$E/\mathrm{eV}$', fontsize = 11)
ax.set_ylabel(r'$\chi/\mathrm{MeVcm}^{2}\mathrm{g}^{-1}$', fontsize = 11)
ax.set_yscale('log')
ax.set_title('Stopping power for protons')
plt.show()
#let's check this
print('huzzah')