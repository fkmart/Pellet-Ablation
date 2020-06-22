import numpy as np 

def k_calcs(x,y,h,BT, func, j, k) :
    if j ==0 :
        k = func(x,y)
    elif j==1:
        k = func(x + BT[1,0]*h, y - h*(k[0]*BT[1,1])) # should be + not - but need to account for negative dr 
    elif j==2:
        k = func(x + BT[2,0]*h, y - h*(k[0]*BT[2,1] + k[1]*BT[2,2]))
    elif j==3:
        k = func(x + BT[3,0]*h, y - h*(k[0]*BT[3,1] + k[1]*BT[3,2] + k[2]*BT[3,3]))
    elif j==4:
        k = func(x + BT[4,0]*h, y - h*(k[0]*BT[4,1] + k[1]*BT[4,2] + k[2]*BT[4,3] + k[3]*BT[4,4]))
    elif j==5:
        k = func(x + BT[5,0]*h, y - h*(k[0]*BT[5,1] + k[1]*BT[5,2] + k[2]*BT[5,3] + k[3]*BT[5,4] + k[4]*BT[5,5])) 
    return k 

def ys_calc(h,k,BT,j):
    if j == 0:
        ys = h*(k[0]*BT[1,1])
    elif j ==1:
        ys = h*(k[0]*BT[2,1] + k[1]*BT[2,2])
    elif j==2:
        ys = h*(k[0]*BT[3,1] + k[1]*BT[3,2] + k[2]*BT[3,3])
    elif j ==3: 
        ys = h*(k[0]*BT[4,1] + k[1]*BT[4,2] + k[2]*BT[4,3] + k[3]*BT[4,4])
    elif j==4:
        ys = h*(k[0]*BT[5,1] + k[1]*BT[5,2] + k[2]*BT[5,3] + k[3]*BT[5,4] + k[4]*BT[5,5])
    return ys