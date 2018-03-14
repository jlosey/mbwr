trange = np.arange(0.7,1.31,0.005)
p = np.array([])
dP = np.array([])
f = open("../data/data1.dat","w")
#for ri in np.arange(0.05,0.4,0.0005):
#    p = np.append(p,mwbrP(ri,0.8))
#    dP = np.append(dP,(ri,0.8))
#print np.diff(p)/0.005
#print dP
for ti in trange:    
    #d1 = newton(dPdV,0.2,args=(ti,),maxiter=400)
    #d2 = newton(dPdV,0.5,args=(ti,),maxiter=400)
    d1 = brentq(dPdV,0.005,0.3,args=(ti,),maxiter=200)
    d2 = brentq(dPdV,0.3,0.9,args=(ti,),maxiter=200)
    #print ti,b,d
    f.write("{}\t{}\t{}\n".format(ti,d1,d2))
f.close()
#plt.plot(rdens,trange):
#plt.show()

