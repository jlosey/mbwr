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

def a1(t,r):
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
    return a

def a2(t,r):
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
    return a

def b1(t,r):
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
    return b

def b2(t,r):
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
    return b

def c1(t,r):
    tsqrt = np.sqrt(t)
    t2 = t**2
    c = []
    c.append(xi[1]*tsqrt/2.0 + xi[2] + 2.0*xi[3]/t + 3.0*xi[4]/t2)
    c.append(xi[6] + 2.0*xi[7]/t + 3.0*xi[8]/t2)
    c.append(xi[10] + 2.0*xi[11]/t)
    c.append(xi[12])
    c.append(2.0*xi[13]/t + 3.0*xi[14]/t2)
    c.append(2.0*xi[15]/t)
    c.append(2.0*xi[16]/t + 3.0*xi[17]/t2)
    c.append(3.0*xi[18]/t2)
    return c

def d1(t,r):
    t2 = t**2
    t3 = t**3
    t4 = t**4
    d = []
    d.append(3.0*xi[19]/t2 + 4.0*xi[20]/t3)
    d.append(3.0*xi[21]/t2 + 5.0*xi[22]/t4)
    d.append(3.0*xi[23]/t2 + 4.0*xi[24]/t3)
    d.append(3.0*xi[25]/t2 + 5.0*xi[26]/t4)
    d.append(3.0*xi[27]/t2 + 4.0*xi[28]/t3)
    d.append(3.0*xi[29]/t2 + 4.0*xi[30]/t3 + 5.0*xi[31]/t4)
    return d

def g1(r):
    gamma = 3.0
    F = math.exp(-1*gamma*r**2)
    g = []
    g.append((1-F)/(2.0*gamma))
    for n in np.arange(2,11,2):
        g.append(-1*(F*r**n-n*g[-1])/(2.0*gamma))
    return g

def a1sum(t,r):
    """Calculate the sum of a1 terms of the MWBR EOS."""
    a = a1(t,r)
    
    total = 0.0
    for i,aa in enumerate(a,1):
        total += aa*r**(i+1)
        
    return total

def b1sum(t,r):
    """Calculate the sum of bi terms in the MWBR EOS"""
    F = math.exp(-3*r**2)
    b = b1(t,r)
    exp = 3.0
    total = 0.0
    for ib in b:
        total += ib*r**exp
        exp += 2
    return total*F

def a2sum(t,r):
    """Calculate the sum of ci terms of the MWBR dP/dT."""
    aii = a2(t,r)
    #ii = 1.0
    total = 0.0
    for ii,ai in enumerate(aii,1):
        total += ai*r**(ii)*(ii+1)
        #ii += 1
    return total

def b2sum(t,r):
    """Calculate the sum of di terms in the MWBR dP/dT"""
    bii = b2(t,r)
    F = math.exp(-3*r**2)
    #ii = 1.0
    total = 0.0
    for ii,bi in enumerate(bii):
        total += bi*r**(2*ii)*(2*ii+1)
        #ii += 1
    return total*F

#def gsum(t,r):
#    """Calculate the density dependent coeffs for the Helmholtz free energy MWBR EOS."""

def mbwrP(dens,tmp):
    """Calculate the pressure for the density and temperature parameters."""
    #X = np.loadtxt("./data/para.dat")
    #print("a1",a1sum(tmp,dens))
    #print("b1",b1sum(tmp,dens))
    press = tmp*dens + a1sum(tmp,dens) + b1sum(tmp,dens)
    return press

def dPdV(dens,tmp):
    """Calculate the pressure-volume derivative."""
    #X = np.loadtxt("./data/para.dat")
    dp = tmp + a2sum(tmp,dens) - 6*dens*b1sum(tmp,dens) + b2sum(tmp,dens)
    return dp

def intrnEng(dens,temp):
    c = c1(temp,dens)
    d = d1(temp,dens)
    g = g1(dens)
    C1 = 0.0
    D1 = 0.0
    for i,ci in enumerate(c,1):
        C1 += (ci*dens**i)/i
    for di,gi in zip(d,g):
        D1 += di*gi
    return C1 + D1

#def cv(dens,tmp):
    
