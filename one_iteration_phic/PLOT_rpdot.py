import numpy as np
import matplotlib.pyplot as plt
import os


l = 500
t = np.linspace(0.0, 1.0,l)
k = [1.0, 5.0, 10.0, 100.0]
k = -np.asarray(k)

title = 'rpdot'

path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))

fig, ax = plt.subplots()
for i in range(0, len(k)):
    rpd = -0.5*k[i]*np.tanh(k[i]*(1.0 - t)) 
    rpd/= np.sqrt(np.log(np.cosh(k[i])))
    rpd /= np.sqrt(np.log(np.cosh(k[i]*(1.0 - t))))
    
    ax.plot(t, rpd, label = r'$k_2 = $' + str(-k[i]))
ax.set_xlabel(r'$\tilde{t}$', fontsize = 13)
ax.set_ylabel(r'$\dot{\tilde{r}}_p$', fontsize = 13, rotation = 0.0)
ax.xaxis.set_label_coords(0.50, -0.04)
ax.yaxis.set_label_coords(-0.05, 0.50)
plt.legend()
fig.tight_layout()
plt.savefig(path + '/' + title + '.png', format = 'png', dpi = 1600)
plt.show()
