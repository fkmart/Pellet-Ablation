import numpy as np
import matplotlib.pyplot as plt 
import mpl_toolkits.mplot3d as m3d 
from gen_var import sig, pel_pot, lp, rc, rp, t_start
import os

direc = os.getcwd() 
load_dir = os.path.join(direc, 'one_iteration_phic', 'analysed_outputs') + os.sep 

fig = plt.figure(figsize = (10.0,8.0))
ax = plt.axes(projection = '3d')

p = 3

pot = format(pel_pot[p], '.1f')
for a in range(1,len(sig)):
    k = format(sig[a], '.2f')
    file = np.loadtxt(load_dir + 'density_pot_test_t'+str(t_start) +'pot'+str(pot) +'sig' + k +'.txt')
    X,Y = np.meshgrid(file[:,1], file[:,0])
    