import numpy as np

def discret(r):
    l = len(r)
    A = np.zeros((l,l))
    # setting up discretisation matrix
    for i in range(0,l):
        A[i,i] = -2.0
        #this could be tidier but not essential to change
    for i in range(1,l):
        A[i,i-1] = 1.0     
    for i in range(0,l-1):
        A[i,i+1] = 1.0  
    return A 