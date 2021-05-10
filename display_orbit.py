import numpy as np
import matplotlib.pyplot as plt

filenames = ['mercpos.txt','venuspos.txt','earthpos.txt','marspos.txt','juppos.txt']
colors = ['salmon','mediumseagreen','magenta','cornflowerblue','limegreen','goldenrod']

x,y,z,xv,yv,zv,t,mass = np.loadtxt('data_6.txt',unpack=True,skiprows=1)
ke= 0.5*mass*(xv**2+yv**2+zv**2)


pe= 1*2.96e-4/np.sqrt((x**2+y**2+z**2))
for k in range(5):
    
    plnt = 'data_'+str(k+1)+'.txt'
    x,y,z,xv,yv,zv,t,mass = np.loadtxt(plnt,unpack=True,skiprows=1)
    ke_arr = 0.5*mass*(xv**2+yv**2+zv**2)
    ke = np.vstack((ke,ke_arr))
    
    pe_arr = 1*2.96e-4/np.sqrt((x**2+y**2+z**2))
    pe = np.vstack((pe,pe_arr))
    plt.scatter(x,y,color=colors[k])
plt.show()
print(pe)
print(ke)
plt.plot(t[10:],pe.sum(axis=0)[10:])
plt.show()
plt.plot(t[10:],ke.sum(axis=0)[10:])
plt.show()
