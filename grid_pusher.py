import numpy as np
import find_nearest as fn 

def pusher(A, r):
    A_out = np.zeros(len(r))
    for i in range(0, len(A)):
        ind = fn.find_nearest(r,A[i,1])
        A_out[ind] += A[i,0]
    return A_out 