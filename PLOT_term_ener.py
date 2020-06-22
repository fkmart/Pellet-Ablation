import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import pel_pot

main_dir = os.getcwd() 
load_dir = os.path.join(main_dir,'one_iteration_phic','analysed_outputs', 't_50') + os.sep 

imp_ener = [] 
for i in range(0, len(pel_pot)):
    load = np.loadtxt(load_dir + 'lin_terminal_energy_pot_' + str(pel_pot[i]) + '_test.txt')
    imp_ener = np.append(imp_ener,load[0,0])

fig,ax = plt.subplots() 
ax.plot(pel_pot, imp_ener)
plt.ylabel(r'$E_{tr}/\mathrm{eV}$', fontsize = 12, rotation = 0)
plt.xlabel(r'$\phi_{\mathrm{pel}}$', fontsize = 12)
ax.xaxis.set_label_coords(0.45,-0.04)
ax.yaxis.set_label_coords(-0.06,0.42)
plt.savefig('min_imp_ener.png', format = 'png', dpi = 1400)
plt.show()