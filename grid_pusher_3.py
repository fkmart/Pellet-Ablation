import numpy as np
import find_nearest as fn 

def pusher(A1,A2, r):
    A_out = np.zeros(len(r))
    for i in range(0, len(A1)):
        ind = fn.find_nearest(r,A1[i])
        A_out[ind] += A2[i]
    return A_out 