import numpy as np
import matplotlib.pyplot as plt 
from gen_var import *
import os 

savedir = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures')
inputdir = './multi_static_outputs/'
term_ener = [] 
filelist = [ f for f in os.listdir(inputdir) if f.startswith("terminal_energy") ]
filelist.sort()
length = len(filelist)
t_plot = t[5::5]
for f in filelist: 
    x = np.loadtxt(inputdir+f)
    term_ener.append(x[0,0])
fig, ax  = plt.subplots()
ax.plot(t_plot, term_ener)

ax.set_xlabel(r'$\tilde{t}$')
ax.set_ylabel('Minimum traversal energy/eV')
ax.xaxis.set_label_coords(0.54, -0.04)
ax.yaxis.set_label_coords(-0.10, 0.45)
plt.savefig(savedir+'/term_ener.png', format = 'png', dpi = 1200)
plt.show()