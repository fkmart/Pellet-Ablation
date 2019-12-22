import numpy as np 

def gcf(bin_arr, x_arr):
    ind = np.argmax(bin_arr)
    max_point = x_arr[ind]
    return max_point