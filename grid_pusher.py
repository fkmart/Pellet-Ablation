import numpy as np
import find_nearest as fn 

def pusher(A, r):
    A_out = np.zeros(len(r))
    for i in range(0, np.shape(A)[1]):
        ind = fn.find_nearest(r,A[1,i])
        A_out[ind] += A[0,i]
    return A_out 