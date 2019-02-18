import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

vo=10
k=np.array([0.5,0.35,0.40,0.21])
#k1=0.5
#k2=0.35
#k3=0.40
#k4=0.21

#Declare the model
def myModel(y,t):
    
    dy=np.zeros(5)
    
    dy[0]= vo -1/k[0]*y[0]
    dy[1]= k[0]*y[0]-k[1]*y[1]
    dy[2]= vo -k[2]*y[2]
    dy[3]= vo-k[0]*y[1]+k[2]*y[2]+y[3]
    dy[4]= k[2]*y[1]-k[3]*y[3]-y[4]
    
    return dy

def solver():

    time=np.linspace(0,20,100)
    yinit=np.zeros(5)
    yinit[4]=1
    y=odeint(myModel,yinit,time)
    
    f, ax = plt.subplots(5) 
    #ax[0,0].plot(time, y[:,0])

    for i in range(0,5):
        ax[i].plot(time, y[:,i])
        ax[i].set_ylabel('y'+str(i))    
        
    plt.xlabel('t')
    plt.suptitle('pyode solver')
    
    plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
    plt.setp(ax[4].get_xticklabels(), visible=True)
    plt.show()