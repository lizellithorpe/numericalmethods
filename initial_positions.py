import numpy as np
import matplotlib.pyplot as plt
from numpy import random
import argparse

import sys
sys.path.append('/Users/lizellithorpe/desktop/research')
from orbital_xyz import *


class MassiveBody:
    def __init__(self,name,mass):
        self.objname = name
        self.mass = mass
        self.x, self.y, self.z = None, None, None
        self.vx, self.vy, self.vz = None, None, None
    
    def assignOrb(self,a,e,inc,capom,om,capm):
        self.a = a
        self.e = e 
        self.inc = inc
        self.capom = capom
        self.om = om 
        self.capm = capm
    
    def assignCart(self):
        self.x,self.y,self.z,self.vx,self.vy,self.vz = orb2xv(np.array([self.a]),np.array([self.e]),
        np.array([self.inc]),np.array([self.capom]),np.array([self.om]),np.array([self.capm]),np.array(self.mass))

class Asteroid:
    def __init__(self):
        self.mass = np.random.random_sample()*(1e15*2.96e-4/1.989e30)
    def assignOrb(self,a,e,inc,capom,om,capm):
        self.a = a
        self.e = e 
        self.inc = inc
        self.capom = capom
        self.om = om 
        self.capm = capm
    
    def assignCart(self):
        self.x,self.y,self.z,self.vx,self.vy,self.vz = orb2xv(np.array([self.a]),np.array([self.e]),
        np.array([self.inc]),np.array([self.capom]),np.array([self.om]),np.array([self.capm]),np.array(self.mass))



txtfile = np.zeros((106,7))
#start with 8 planets and 10 test particle asteroid belt

mercury = MassiveBody('mercury',2.96e-4*1.66e-7)
mercury.assignOrb(0.387,0.205,7.004*np.pi/180,48.33*np.pi/180,77.45*np.pi/180,np.random.random_sample()*2*np.pi)
mercury.assignCart()

txtfile[:][0] = np.array([mercury.x,mercury.y,mercury.z,mercury.vx,mercury.vy,mercury.vz,mercury.mass])


venus = MassiveBody('venus',2.96e-4*2.45e-6)
venus.assignOrb(0.723,0.006,3.394*np.pi/180,76.68*np.pi/180,131.53*np.pi/180,np.random.random_sample()*2*np.pi)
venus.assignCart()

txtfile[:][1] = np.array([venus.x,venus.y,venus.z,venus.vx,venus.vy,venus.vz,venus.mass])
earth = MassiveBody('earth',3.003e-6)
earth.assignOrb(1,0.016,0.00005*np.pi/180,-1.26*np.pi/180,102.94*np.pi/180,np.random.random_sample()*2*np.pi)
earth.assignCart()
txtfile[:][2] = np.array([earth.x,earth.y,earth.z,earth.vx,earth.vy,earth.vz,earth.mass])
mars = MassiveBody('mars',3.23e-7)
mars.assignOrb(1.52,0.093,1.85*np.pi/180,49.57*np.pi/180,336.04*np.pi/180,np.random.random_sample()*2*np.pi)
mars.assignCart()
txtfile[:][3] = np.array([mars.x,mars.y,mars.z,mars.vx,mars.vy,mars.vz,mars.mass])
jupiter = MassiveBody('jupiter',2.96e-4*0.0009543)
jupiter.assignOrb(5.203,0.048,1.305*np.pi/180,100.556*np.pi/180,14.75*np.pi/180,np.random.random_sample()*2*np.pi)
jupiter.assignCart()


txtfile[:][4] = np.array([jupiter.x,jupiter.y,jupiter.z,jupiter.vx,jupiter.vy,jupiter.vz,jupiter.mass])
print(txtfile)
#make some random asteroids
objs = [Asteroid() for j in range(100)]

for k,obj in enumerate(objs):
    obj.assignOrb(2.06+np.random.random_sample()*1.21,np.random.random_sample()*0.33,np.random.random_sample()*20*np.pi/180,
    np.random.random_sample()*2*np.pi,np.random.random_sample()*2*np.pi,np.random.random_sample()*2*np.pi)
    obj.assignCart()
    txtfile[:][k+6] = np.array([obj.x,obj.y,obj.z,obj.vx,obj.vy,obj.vz,obj.mass])
txtfile[:][5] = np.array([0,0,0,0,0,0,2.96e-4])
print(txtfile)
np.savetxt('initial_pos.txt',txtfile)
#write initial positions, velocities into txt file
