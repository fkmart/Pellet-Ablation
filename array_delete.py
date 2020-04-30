import numpy as np

def array_delete_less(A,val):
    arr = []
    for i in range(0, len(A)):
        if A[i] < val:
            arr = np.append(arr,i)
        else:
            pass
    return arr

def array_delete_more(A,val):
    arr = []
    for i in range(0, len(A)):
        if A[i] > val:
            arr = np.append(arr,i)
        else:
            pass
    return arr