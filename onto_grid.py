import numpy as np 

"""This function is designed to take the stopping points of the electrons
and shift those electrons from their VERY fine resolution positions determined by
a NOT finite difference RK4 algorithm to their required finite difference positions
on an appropriately resolved grid for SOR.
""" 

def onto_grid(stop_points, dens, fd_grid, grid_dens):  
    def find_nearest(arr,val):
        ind = (np.abs(arr[:] - val)).argmin()
        return ind
    for i in range(0, len(stop_points)):
        ind = find_nearest(fd_grid, stop_points[i])
        grid_dens[ind] += dens[i]
    return grid_dens

points = np.linspace(0.0, 1.0, num =200)
dens = np.linspace(0.0, 0.5, num = 200)
fd_grid = np.arange(0.0, 1.0, 0.02)
grid_dens = np.zeros(len(fd_grid))

grid_dens = onto_grid(points,dens, fd_grid, grid_dens)
print(grid_dens)