# -*- coding: utf-8 -*-
'''Python code developed running Canopy (Free Academic License)
https://www.enthought.com/products/canopy/'''
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
#np.set_printoptions(threshold=np.inf)

def vectorfielda(w, t, par):
    
    """
    *Problem 8-05 Part (d)
    {Elements of Chemical Reaction Engineering, Fogler (5th Ed.)}
    
    Arguments:
        w :  vector of the state variables:
                  w = [Ca,Cb]
        W :  catalyst weight
        par :  vector of the parameters:
                  par = [k1,k2]
                  
                  A--k1--> B --k2--> C
    """
 
    Ca, Cb = w
    k1,k2 = par
    
    #Drinks consumed at a continuous rate 1st hour (change of conc. due to incoming EtOH is 2 g L-1 h-1)
    # Create f = (dCa/dt,dCb/dt)
    f = [-k1*Ca, k1*Ca-k2]

    return f
    
def vectorfieldb(w, t, par):
    
    Ca, Cb = w
    k1,k2 = par
    
    #EtOH absorbed and cleared remaining hours
    # Create f = (dCa/dt,dCb/dt)
    f = [-k1*Ca, k1*Ca-k2]

    return f
    
    
def solver():
    # Use ODEINT to solve the differential equations defined by the vector field
    #from scipy.integrate import odeint
    
    k1=10
    k2=0.192

    
    # ODE solver parameters
    abserr = 1.0e-8
    relerr = 1.0e-6
   
    stopxa = 0.5
    numpointsa = 80
    stopxb = 10
    numpointsb = 80
    
    # Create the x-axis samples for the output of the ODE solver.
    ta = [stopxa * float(i) / (numpointsa - 1) for i in range(numpointsa)]
    tb = [stopxa + (stopxb * float(i) / (numpointsb - 1)) for i in range(numpointsb)]

    # Initial conditions before 1st hour
    Ca0 = 1
    Cb0 = 0
    
    # Pack up the parameters and initial conditions:
    par = [k1,k2]
    w0 = [Ca0,Cb0]
    
    # Call the ODE solver.
    wsola = odeint(vectorfielda, w0, ta, args=(par,), atol=abserr, rtol=relerr)
    
    
    # Initial conditions after 1st hour
    Ca02 = 1
    Cb02 = wsola[79,1]
    w02 = [Ca02,Cb02]
    
    wsolb = odeint(vectorfieldb, w02, tb, args=(par,), atol=abserr, rtol=relerr)

 
    with open('BAC805ad.dat', 'w') as f:
        # Print & save the solution.
        for t1, w1 in zip(ta, wsola):
            print >> f, t1, w1[0], w1[1]
    
    with open('BAC805bd.dat', 'w') as f:
        # Print & save the solution.
        for t2, w2 in zip(tb, wsolb):
            print >> f, t2, w2[0], w2[1]
        
    # Plot the solution that was generated
    
    from numpy import loadtxt
    from pylab import savefig
    from matplotlib.font_manager import FontProperties
    
    t1, Ca1, Cb1 = loadtxt('BAC805ad.dat', unpack=True)    
    t2, Ca2, Cb2 = loadtxt('BAC805bd.dat', unpack=True)

    
    lw = 1
    
    # Two subplots, unpack the axes array immediately
    f, (ax1) = plt.subplots(1, 1, sharex=True)
    #grid(True)
    ax1.plot(t1, Ca1, 'r',linestyle='--', linewidth=2)
    ax1.plot(t1, Cb1, 'b',linestyle='--', linewidth=lw)
    ax1.plot(t2, Ca2, 'r', linewidth=2)
    ax1.plot(t2, Cb2, 'b', linewidth=lw)
    ax1.legend((r'$C_A$', r'$C_B$',r'$C_{A,after 1/2h}$', r'$C_{B,after 1/2h}$'), prop=FontProperties(size=16))
    #ax1.set_title('')
    ax1.set_xlabel('$t$')
    ax1.set_ylim(0)
    
    
    f.subplots_adjust(hspace=0,wspace=0.19,left=0.06,right=0.96,top=0.94,bottom=.10)
    
    savefig('BAC805d.png', dpi=100) 
