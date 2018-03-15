import numpy as np
import math

xi = np.array([
    0.8623085097507421,
    2.976218765822098,
    -8.402230115796038,
    0.1054136629203555,
    -0.8564583828174598,
    1.582759470107601,
    0.7639421948305453,
    1.753173414312048,
    2.798291772190376e+03,
    -4.8394220260857657e-02,
    0.9963265197721935,
    -3.698000291272493e+01,
    2.084012299434647e+01,
    8.305402124717285e+01,
    -9.574799715203068e+02,
    -1.477746229234994e+02,
    6.39860785271505e+01,
    1.603993673294834e+01,
    6.805916615864377e+01,
    -2.791293578795945e+03,
    -6.245128304568454,
    -8.116836104958410e+03,
    1.488735559561229e+01,
    -1.059346754655084e+04,
    -1.131607632802822e+02,
    -8.867771540418822e+03,
    -3.986982844450543e+01,
    -4.689270299917261e+03,
    2.593535277438717e+02,
    -2.694523589434903e+03,
    -7.218487631550215e+02,
    1.721802063863269e+02])

def asum(t,r):
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

def bsum(t,r):
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

def mbwrP(dens,tmp):
    """Calculate the pressure for the density and temperature parameters."""
    #X = np.loadtxt("./data/para.dat")
    press = tmp*dens + asum(tmp,dens) + bsum(tmp,dens)
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

