import numpy as np 

def find_nearest(A,val):
    ind = (np.abs(A[:] - val)).argmax()
    return ind 
    