import numpy as np
import math

def asum(xi,t,r):
    """Calculate the sum of ai terms of the MWBR EOS."""
    tsqrt = np.sqrt(t)
    t2 = t**2
    a = []
    a.append(xi[0]*t + xi[1]*tsqrt + xi[2] + xi[3]/t + xi[4]/t2)
    a.append(xi[5]*t + xi[6] + xi[7]/t + xi[8]/t2)
    a.append(xi[9]*t + xi[10] + xi[11]/t)
    a.append(xi[12])
    a.append(xi[13]/t + xi[14]/t2)
    a.append(xi[15]/t)
    a.append(xi[16]/t + xi[17]/t2)
    a.append(xi[18]/t2)
    exp = 1.0
    total = 0.0
    for ai in a:
        total += ai*r**(exp+1)
        exp += 1
    return total

def bsum(xi,t,r):
    """Calculate the sum of bi terms in the MWBR EOS"""
    F = math.exp(-3*r**2)
    t2 = t**2
    t3 = t**3
    t4 = t**4
    b = []
    b.append(xi[19]/t2 + xi[20]/t3)
    b.append(xi[21]/t2 + xi[22]/t4)
    b.append(xi[23]/t2 + xi[24]/t3)
    b.append(xi[25]/t2 + xi[26]/t4)
    b.append(xi[27]/t2 + xi[28]/t3)
    b.append(xi[29]/t2 + xi[30]/t3 + xi[31]/t4)
    exp = 3.0
    total = 0.0
    for bi in b:
        total += bi*r**exp
        exp += 2
    return total*F

def mwbrP(dens,tmp):
    X = np.loadtxt("./data/para2.dat")
    press = tmp*dens + asum(X,tmp,dens) + bsum(X,tmp,dens)
    return press

def csum(xi,t,r):
    """Calculate the sum of ci terms of the MWBR dP/dT."""
    tsqrt = np.sqrt(t)
    t2 = t**2
    a = []
    a.append(xi[0]*t + xi[1]*tsqrt + xi[2] + xi[3]/t + xi[4]/t2)
    a.append(xi[5]*t + xi[6] + xi[7]/t + xi[8]/t2)
    a.append(xi[9]*t + xi[10] + xi[11]/t)
    a.append(xi[12])
    a.append(xi[13]/t + xi[14]/t2)
    a.append(xi[15]/t)
    a.append(xi[16]/t + xi[17]/t2)
    a.append(xi[18]/t2)
    ii = 1.0
    total = 0.0
    for ai in a:
        total += ai*r**(ii)*(ii+1)
        ii += 1
    return total

def dsum(xi,t,r):
    """Calculate the sum of di terms in the MWBR dP/dT"""
    F = math.exp(-3*r**2)
    t2 = t**2
    t3 = t**3
    t4 = t**4
    b = []
    b.append(xi[19]/t2 + xi[20]/t3)
    b.append(xi[21]/t2 + xi[22]/t4)
    b.append(xi[23]/t2 + xi[24]/t3)
    b.append(xi[25]/t2 + xi[26]/t4)
    b.append(xi[27]/t2 + xi[28]/t3)
    b.append(xi[29]/t2 + xi[30]/t3 + xi[31]/t4)
    ii = 1.0
    total = 0.0
    for bi in b:
        total += bi*r**(2*ii)*(2*ii+1)
        ii += 1
    return total*F

def dPdV(dens,tmp):
    X = np.loadtxt("./data/para.dat")
    dp = tmp + csum(X,tmp,dens) - 6*dens*bsum(X,tmp,dens) + dsum(X,tmp,dens)
    return dp

