#!/usr/bin/env python

import math
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from scipy import signal

def iterateU(Uarray, Harray, index):
    global g
    global n
    global t
    global x

    a = 0
    
    if index==0:
        a  = Uarray[index] - g*(t/x)*(Harray[index] - 0)
    else:
        a = Uarray[index] - g*(t/x)*(Harray[index]-Harray[index-1])
    
    return a

def iterateH(Uarraynext, Harray, Barray, index):
    global n
    global t
    global x

    a = 0

    if index==n-1:
        a = Harray[index] - Barray[index]*(t/x)*(0 - Uarraynext[index])
    else:
        a = Harray[index] - Barray[index]*(t/x)*(Uarraynext[index+1]-Uarraynext[index])

    return a

def smooth(arrayh):
    global n

    a = np.zeros(n)

    for i in range(n-1):
        a[i] += (arrayh[i+1]-arrayh[i])*0.2*0.5
        a[i+1] +=  (arrayh[i+1]-arrayh[i])*0.2*0.5*(-1)

    for i in range(n):
        arrayh[i] += a[i]

    return arrayh

if __name__ == "__main__":
    n = 1000     #number of bins
    x = 111.0     #distance between bins
    t = 10    #time step between iterations
    g = 9.81     #strength of gravity
    d = 1000.0     #depth of water
    timesteps = 100     #number of iterations
    passes = 2    #number of smoothing passes
    tmax = x/(math.sqrt(g*d))     #maximum time step allowed

    print("maximum t: {}".format(tmax))

    if (t>tmax):
        t=tmax
        print("new t: {}".format(t))
    
    h = np.zeros(n)
    u = np.zeros(n)

    for i in range(n):
        if i<(n/2+10) and i>(n/2-10):
            h[i] = 1
        if (i<=(n/2-10) and i>(n/2-20)) or (i>=(n/2+10) and i<(n/2+20)):
            h[i] = 0

    bathymetry = np.zeros(n)
    for i in range(n):
        if i>n*0.75:
            bathymetry[i] = d
        else:
            bathymetry[i] = d
        
    finalh = np.zeros([timesteps+1,n])
    finalh[0] = h

    finalu = np.zeros([timesteps+1,n])
    finalu[0] = u

    a = np.zeros(n)
    b = np.zeros(n)
    e = np.zeros(n)
    
    for i in range(timesteps):
        for j in range(n):
            a[j] = iterateU(finalu[i], finalh[i], j)
            
        finalu[i+1] = a

        for k in range(passes):
            finalu[i+1] = smooth(finalu[i+1])
                
        for j in range(n):
            a[j] = iterateH(finalu[i+1], finalh[i], bathymetry, j)
            
        finalh[i+1] = a

    xaxis = np.zeros(n)
    for i in range(n):
        xaxis[i]=i*x
               
    fig, ax = plt.subplots()
    ln, = plt.plot([], [], 'b-', animated=True)
    
    def init():
        ax.set_xlim(0,(n-1)*x)
        ax.set_ylim(-1,1)
        return ln,
        
    def update(frame):
        ln.set_data(xaxis, finalh[frame])
        return ln,
        
    ani = FuncAnimation(fig, update, frames=timesteps+1, init_func=init, interval=100, blit=True)
    
    plt.show()

    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    
    ani.save('waveanimation.mp4', writer='ffmpeg')
