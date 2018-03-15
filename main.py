from mbwr import eos 
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

trange = np.arange(0.7,1.31,0.005)
p = np.array([])
dP = np.array([])
#f = open("../data/data1.dat","w")
#for ri in np.arange(0.05,0.4,0.0005):
#    p = np.append(p,mwbrP(ri,0.8))
#    dP = np.append(dP,(ri,0.8))
#print np.diff(p)/0.005
for ti in trange:    
    #d1 = newton(dPdV,0.2,args=(ti,),maxiter=400)
    #d2 = newton(dPdV,0.5,args=(ti,),maxiter=400)
    d1 = brentq(eos.dPdV,0.005,0.3,args=(ti,),maxiter=200)
    d2 = brentq(eos.dPdV,0.3,0.9,args=(ti,),maxiter=200)
    print ti,d1,d2
    #f.write("{}\t{}\t{}\n".format(ti,d1,d2))
#f.close()
#plt.plot(d1,trange,"-r",d2,trange,"-k")
#plt.show()

