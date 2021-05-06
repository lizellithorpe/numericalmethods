import numpy as np
import matplotlib.pyplot as plt

filenames = ['mercpos.txt','venuspos.txt','earthpos.txt','marspos.txt','juppos.txt']
colors = ['salmon','mediumseagreen','magenta','cornflowerblue','limegreen']
for k,plnt in enumerate(filenames):
    print(plnt)
    x,y,z,t,mass = np.loadtxt(plnt,unpack=True)
    plt.scatter(x,y,color=colors[k])
plt.show()