import numpy as np 
"""
a = np.arange(0,10,1)
b = np.zeros(len(a))
b[0], b[1] = 1.0, 0.5
b[2::2] = 0.25
b[3::2] = 0.125
c = 0 
z = []
for i in range(0, len(a)):
    c +=1
    if c%8 == 0:
        x = b[c-8:c]
        y = np.sort(x)
        for j in range(0,len(x)-1):
            if y[j] == y[j+1]:
                z.append(y[j])
        if len(z) >0:
            z = np.unique(z)
            ind1 = np.where(x == np.amax(z))[0]
            ind2 = np.where(x == np.amin(z))[0]
            l1 = len(ind1)
            l2 = len(ind2)
            if l1 > l2:
                ind1 = ind1[:l2]
            elif l2 > l1:
                ind2 = ind2[:l1]
            else:
                pass 
            for j in range(0, len(ind1)):
                if (ind1[j] == ind2[j] -1) or (ind1[j] == ind2[j] +1):
                    d = np.amin(z)
                else:
                    break                
        else:
            pass
"""
def frigg_it_logic(b,c):
    z = []
    #x = b[c-8:c]
    x = b[-c:]
    y = np.sort(x)
    for j in range(0,len(x)-1):
        if y[j] == y[j+1]:
            z.append(y[j])
    if len(z) >0:
        z = np.unique(z)
        ind1 = np.where(x == np.amax(z))[0]
        ind2 = np.where(x == np.amin(z))[0]
        l1 = len(ind1)
        l2 = len(ind2)
        if l1 > l2:
            ind1 = ind1[:l2]
        elif l2 > l1:
            ind2 = ind2[:l1]
        else:
            pass 
        for j in range(0, len(ind1)):
            if (ind1[j] == ind2[j] -1) or (ind1[j] == ind2[j] +1):
                d = 'true'
            else:
                d = 'false'
                break                
        else:
            pass
    step = np.amin(z)
    return d, step