#!/usr/bin/env python

import math
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from scipy import signal
from netCDF4 import Dataset

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

def smooth(arrayhu):
    global n

    a = np.zeros(n)

    for i in range(n-1):
        a[i] += (arrayhu[i+1]-arrayhu[i])*0.15*0.05
        a[i+1] +=  (arrayhu[i+1]-arrayhu[i])*0.15*0.05*(-1)

    for i in range(n):
        arrayhu[i] += a[i]

    return arrayhu

def gaus(mean, sd, volume, where):
    
    return volume*(1/(sd*math.sqrt(2*math.pi)))*math.exp(-((where-mean)**2)/(2*sd**2))
    
    
    
if __name__ == "__main__":
    n = 1000     #number of bins
    x = 100.0     #distance between bins
    t = 2.0    #time step between iterations
    g = 9.81     #strength of gravity
    d = 10.0     #depth of water
    timesteps = 500     #number of iterations
    passes = 0    #number of smoothing passes
    tmax = x/(math.sqrt(g*d))     #maximum time step allowed

    print("maximum t: {}".format(tmax))

    if (t>tmax):
        t=tmax
        print("new t: {}".format(t))
    
    h = np.zeros(n)
    u = np.zeros(n)

    vol = 10*0.1*x
    mean = float(x*(n/2+1))
    sd = 10*x/2
    for i in range(n):
            h[i] = gaus(mean, sd, vol, i*x)
            print(h[i],mean,sd,vol,i*x)

    bathymetry = np.zeros(n)
    for i in range(n):
        if i>n*0.6:
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

        for k in range(passes):
            finalh[i+1] = smooth(finalh[i+1])
            
    xaxis = np.zeros(n)
    for i in range(n):
        xaxis[i]=i*x
               
    fig, ax = plt.subplots()
    ln, = plt.plot([], [], 'b-', animated=True)
    
    def init():
        ax.set_xlim((n-1)*x*3/8,(n-1)*x*5/8)
        ax.set_ylim(-0.02,0.1)
        return ln,
        
    def update(frame):
        ln.set_data(xaxis, finalh[frame])
        return ln,
        
    ani = FuncAnimation(fig, update, frames=timesteps+1, init_func=init, interval=10, blit=True)

    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title='MySimulation', artist='Matplotlib', comment='Animation')
    writer = FFMpegWriter(fps=15, metadata=metadata, bitrate=1000)

    ani.save("gaussianinitialcondition.mp4", writer=writer)
    
    plt.show()
