import numpy as np

a = np.linspace(0.0,  5.0, endpoint = 'true', num = 6)

for i in range(0, len(a)):
    c = 0
    while (c < 4):
        a[i] += 1.0
        if (a[i] > 7.0):
            break
        else:
            pass
        c +=1
    print(str(i)+ 'th iteration is complete')
