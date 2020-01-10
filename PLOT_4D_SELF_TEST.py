import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm

x = np.linspace(-5,6,22)
y = np.linspace(-3,4,14)

X,Y = np.meshgrid(x,y)
z = np.ones((y.size,x.size,3))
z[:,:,1] *= 2
z[:,:,2] *= 3 

A = np.zeros(z.shape)

fig = plt.figure()
ax = fig.gca(projection='3d')
c = 2
for i in range(0,3):
    arr = np.linspace(-c,c+1,x.size)
    array, byp = np.meshgrid(arr,y)
    A[:,:,i] = array
    c+=1
vmi = np.min(A)
vma = np.max(A)

norm = matplotlib.colors.Normalize(vmin = vmi,vmax = vma)
xs = X.shape
ys = Y.shape
zs = z.shape
for i in range(0, 3):
    ax.plot_surface(X,Y,z[:,:,i],facecolors = cm.viridis(norm(A[:,:,i])))
m = cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
m.set_array([])
plt.colorbar(m)
plt.show()