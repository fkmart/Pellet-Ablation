import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.colors as colors
from electron import M_fac, RME
from gen_var import r0cgs, solid_dens,I, zovera, eps

E = np.arange(100,20000, 100)

nondim = (r0cgs/RME)*solid_dens
c1 = 0.153536 
D = np.logspace(-2,1,num = 100, endpoint = 'true')

X,Y = np.meshgrid(E,D)
T = X/(RME*M_fac)
B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula            
S = Y*c1*zovera*B/(2.0*T) #final stopping power form

fig, ax = plt.subplots() 

p = ax.pcolor(E,D,S,norm = colors.LogNorm(vmin = S.min(), vmax = S.max()),cmap = 'gnuplot')

#im = plt.imshow(S, interpolation = 'bilinear', cmap = cm.RdYlGn,
# origin = 'lower', vmax = S.max(), vmin = S.min())
fig.colorbar(p, ax=ax, label = r'$S/\mathrm{MeVcm}^{-1}$')
plt.xlabel(r'$E/\mathrm{eV}$')
plt.ylabel(r'$\rho / \mathrm{gcm}^{-3}$')
plt.yscale('log')
plt.savefig('bethe_colour_plot.png', format = 'png', dpi = 1400)
plt.show()
