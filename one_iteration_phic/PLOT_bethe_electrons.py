

import numpy as np
import os
from gen_var import *
from electron import RME, M_fac
import matplotlib.pyplot as plt
E = np.arange(50,20000,100)

title = 'bethe electrons'

path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))
T = E/(RME*M_fac)

den = 1.0
nondim = (r0cgs/RME)*solid_dens #non-dimensionalisation factor, in cgs to work with formula,this eliminates units of c1
#nondim = 1.0 #THIS IS TO TEST WHETHER NONDIMENSIONALITY ACTUALLY MAKES A DIFFERENCE
c1= 0.153536 #factor that must be an amalgamation of many constants
B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula            
dxde_gas = nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
dxde_real = dxde_gas *(r0cgs/(RME))**(-1)

fig, ax = plt.subplots()
ax.plot(E,dxde_real)
ax.set_xlabel(r'$E/\mathrm{eV}$', fontsize = 12)
ax.set_ylabel('$S$/ \n $\mathrm{MeVcm}^{-1}$', fontsize = 12, rotation = 0 )
ax.yaxis.set_label_coords(-0.09, 0.55)
ax.set_yscale('log')
fig.tight_layout()
#ax.set_title('Stopping power for electrons')
plt.savefig(path + '/' + title + '.png', format = 'png', dpi = 1600)
plt.show()
#let's check this
print('huzzah')