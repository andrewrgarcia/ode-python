# -*- coding: utf-8 -*-
'''Python code developed running Canopy (Free Academic License)
https://www.enthought.com/products/canopy/'''
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
#np.set_printoptions(threshold=np.inf)

def vectorfield(w, t, par):
    
    """
    *Problem 8-05 Part (a)
    {Elements of Chemical Reaction Engineering, Fogler (5th Ed.)}
    
    Arguments:
        w :  vector of the state variables:
                  w = [Ca,Cb]
        t :  time
        par :  vector of the parameters:
                  par = [k1,k2]
    """

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

    stopx = 11
    numpoints = 80

    # Create the x-axis samples for the output of the ODE solver.
    t = [stopx * float(i) / (numpoints - 1) for i in range(numpoints)]

    # Initial conditions before 1st hour
    Ca0 = 2
    Cb0 = 0
    
    # Pack up the parameters and initial conditions:
    par = [k1,k2]
    w0 = [Ca0,Cb0]
    
    # Call the ODE solver.
    wsol = odeint(vectorfield, w0, t, args=(par,), atol=abserr, rtol=relerr)


 
    with open('BAC805x.dat', 'w') as f:
        # Print & save the solution.
        for t1, w1 in zip(t, wsol):
            print >> f, t1, w1[0], w1[1]
    
        
    # Plot the solution that was generated
    
    from numpy import loadtxt
    from pylab import savefig
    from matplotlib.font_manager import FontProperties
    
    t, Ca, Cb = loadtxt('BAC805x.dat', unpack=True)    

    lw = 1
    
    # Two subplots, unpack the axes array immediately
    f, (ax1) = plt.subplots(1, 1, sharex=True)
    #grid(True)
    ax1.plot(t, Ca, 'r', linewidth=2)
    ax1.plot(t, Cb, 'b', linewidth=lw)
    ax1.legend((r'$C_A$', r'$C_B$'), prop=FontProperties(size=16))
    #ax1.set_title('')
    ax1.set_xlabel('$t$')
    ax1.set_ylim(0)
    
    
    f.subplots_adjust(hspace=0,wspace=0.19,left=0.06,right=0.96,top=0.94,bottom=.10)
    
    savefig('BAC805a.png', dpi=100) 
