import numpy as np
import matplotlib.pyplot as plt 
import os 

cwd = os.getcwd()
load_path = os.path.join(cwd,'many_iteration_TEST','analysed_outputs') + os.sep
many_start = 50
fig, ax = plt.subplots(figsize = (10,8))
ax2 = ax.twinx()
lines = []
for i in range(many_start, many_start + 40,10):
    file = np.loadtxt(load_path + 'density_pot_test_t' + str(i) + '.txt')
    ax.scatter(file[:,0], file[:,1],  marker = 'x')
    
for i in range(many_start, many_start + 40,10):    
    file2 = np.loadtxt(load_path + 'potential_update_test_t' + str(i) + '.txt')
    time = format(i/500,'.2f')
    lin = ax2.plot(file2[:,1], file2[:,0],label = r'$\tilde{t} = $' + str(time))
    

ax.set_ylabel('EEDF Bins', rotation = 0)
ax.yaxis.set_label_coords(-0.05, 0.45)
ax.set_xlabel(r'$\tilde{r}$')
ax2.set_ylabel(r'$\phi$' + '/V', rotation = 0)
plt.text(0.5,-500.0, r'$\leftarrow$ to pellet')
plt.text(8.0,-500.0, r'$\rightarrow$ to plasma')
plt.legend(loc = 'center left')
plt.savefig('birthing_gaussians_bins.png', format = 'png', dpi = 1200)
plt.show()