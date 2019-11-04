import numpy as np 
import matplotlib.pyplot as plt  
import math as mt 

savedir = '/home/kyle/Documents/thesis/Local_figures'
title = 'dadt_D2'

D2 = -1.0 
k = -1.0*np.asarray([1.0, 5.0, 10.0, 100.0])
t = np.linspace(0, 1.0, 1000)

dadt = np.zeros((len(t), len(k)))
for i in range(0,len(k)):
    for j in range(0, len(t)):
        dadt[j,i] = mt.tanh(k[i]*(1.0 - t[j]))

fig, ax  = plt.subplots() 
lines = []

for i in range(0, len(k)):
    line = ax.plot(t, dadt[:,i], label = r'$k_2 = $' + str(-k[i]))
    lines += line

line_flat = ax.plot(t, -np.ones(len(t)), linestyle = '--', color = 'black', label = r'$D^2$')
labels = [line.get_label() for line in lines]
labels.append(line_flat[0].get_label())

ax.set_xlabel(r'$\tilde{t}$', fontsize = 14)
ax.set_ylabel(r'$\frac{\mathrm{d}\tilde{A}}{\mathrm{d}\tilde{t}}$', rotation  = 0, fontsize = 14)

ax.yaxis.set_label_coords(-0.04, 0.45)
ax.xaxis.set_label_coords(0.52, -0.05)
plt.legend(labels)
plt.savefig(savedir+'/'+title+'.png', format = 'png', dpi = 1200)
plt.show()