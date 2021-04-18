import numpy as np
import matplotlib.pyplot as plt

filenames = ['mercpos.txt','venuspos.txt','earthpos.txt','marspos.txt','juppos.txt']

for plnt in filenames:
    print(plnt)
    x,y,z,t,mass = np.loadtxt(plnt,unpack=True)
    plt.scatter(x,y)
plt.show()