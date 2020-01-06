import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

x = np.linspace(0.0, 5.0, 100)
y = np.linspace(0.0, 2.0, 50)
X,Y = np.meshgrid(x,y)
Z = np.cos(X) + np.cos(Y)

C = np.cos(np.pi*Y)


color_dimension = C # change to desired fourth dimension
minn, maxx = color_dimension.min(), color_dimension.max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='plasma')
m.set_array([])
fcolors = m.to_rgba(color_dimension)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X,Y,Z, rstride=1, cstride=1, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
fig.colorbar(m)
plt.show()