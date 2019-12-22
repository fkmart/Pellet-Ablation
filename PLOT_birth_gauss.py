import numpy as np
import matplotlib.pyplot as plt 
from gen_var import rp, rc
from electron import RME, M_fac
import os 

cwd = os.getcwd()
load_path = os.path.join(cwd,'many_iteration_TEST','analysed_outputs') + os.sep
many_start = 50
fig, ax = plt.subplots(figsize = (10.0,8.0))
lines = []
colors = ['navy', 'forestgreen', 'gold', 'firebrick']
c = 0
lines += ax.plot([],[], label = "Pre-flow", color = 'white')
for i in range(many_start, many_start + 40,10):
    file = np.loadtxt(load_path + 'density_pot_test_t' + str(i) + '.txt')
    file1 = np.loadtxt(load_path + 'potential_update_test_t' + str(i) + '.txt')
    time = format(i/500,'.2f')
    lines += ax.plot(file1[:,1], file1[:,0], linestyle = '--', color = colors[c],
        label = r'$\tilde{t} = $' + str(time) )
    c +=1

c = 0
lines +=ax.plot([],[], label = "Post-flow", color = 'white')
for i in range(many_start, many_start + 40,10):    
    file2 = np.loadtxt(load_path + 'advected_potential_t' + str(i) + '.txt')
    time = format(i/500,'.2f')
    time2 = format((i+10)/500 ,'.2f')
    lines += ax.plot(file2[:,1], file2[:,0],label = r'$\tilde{t} = $' + str(time) + r'$\rightarrow \tilde{t} = $' + str(time2), color = colors[c])
    c+=1
    
labels = [l.get_label() for l in lines]
ax.set_ylabel('Potential/V', rotation = 0)
ax.yaxis.set_label_coords(-0.05, 0.45)
ax.set_xlabel(r'$\tilde{r}$')
#ax2.set_ylabel(r'$\phi$' + '/V', rotation = 0)
plt.text(0.5,-500.0, r'$\leftarrow$ to pellet')
plt.text(8.0,-500.0, r'$\rightarrow$ to plasma')
plt.legend(labels,loc = 'center left')
#leg = plt.legend()
#for lines,text in zip(leg.get_lines(), leg.get_texts()):
    #if text.get_color() == 'white':
     #   text.set_weight('bold')
plt.grid('true')
plt.savefig('birthing_gaussians_pot_flow.png', format = 'png', dpi = 1400)
plt.show()