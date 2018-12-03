#!/usr/bin/env python

import math
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from scipy import signal
from netCDF4 import Dataset

#iterates velocity array
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

#iterates height array
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

#smoothes height or velocity array
def smooth(arrayhu):
    global n

    a = np.zeros(n)

    for i in range(n-1):
        a[i] += (arrayhu[i+1]-arrayhu[i])*0.15*0.5
        a[i+1] +=  (arrayhu[i+1]-arrayhu[i])*0.15*0.5*(-1)

    for i in range(n):
        arrayhu[i] += a[i]

    return arrayhu

if __name__ == "__main__":

    #initialize variables
    n = 1000     #number of bins
    x = 5585.0     #distance between bins
    t = 100    #time step between iterations
    g = 9.81     #strength of gravity
    d = 1000.0     #depth of water
    timesteps = 800     #number of iterations
    passes = 0    #number of smoothing passes
    tmax = x/(math.sqrt(g*d))     #maximum time step allowed

    print("maximum t: {}".format(tmax))

    if (t>tmax):
        t=tmax
        print("new t: {}".format(t))

    #creates initial height and velocity arrays
    h = np.zeros(n)
    u = np.zeros(n)

    waveheight = 1.0
    wavelength = x*100
    wavenumber = float(2*math.pi/wavelength)
    wavevelocity = math.sqrt(g*d)
    
    for i in range(n):
        if i>=200:
            h[i] = 0
        else:
            h[i] = 2*waveheight*math.sin(-wavenumber*i*x)

    for i in range(n):
        if i<n*0.5:
            u[i] = 0 #math.sqrt(g*d)

    d2 = d/2
    bathymetry = np.zeros(n)
    for i in range(n):
        if i<n*0.5:
            bathymetry[i] = d
        else:
            bathymetry[i] = d2
        
    finalh = np.zeros([timesteps+1,n])
    finalh[0] = h

    finalu = np.zeros([timesteps+1,n])
    finalu[0] = u

    #iterates over all of the bins to create height and velocity arrays for each timestep
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

    #creating x-axis arrays
    xaxis = np.zeros(n)
    xaxis2 = np.zeros(n)
    for i in range(n):
        xaxis[i]=i*x
        xaxis2[i]=i*x

    #calculates reflected and transmitted factors
    reflectedFactor = -(math.sqrt(d) - math.sqrt(d2))/(math.sqrt(d)+math.sqrt(d2))
    transmittedFactor = 1 + reflectedFactor
    
    #animates arrays
    fig, ax = plt.subplots()
    ln, = plt.plot([], [], 'b-', linewidth=4, animated=True)
    ln2, = plt.plot([], [], 'r-', animated=True)
    
    def init():
        ax.set_xlim(0,(n-1)*x)
        ax.set_ylim(-1.5,1.5)
        return ln, ln2, 
        
    def update(frame):
        ln.set_data(xaxis, finalh[frame])
        for i in range(n):
            if i<n*0.5:
                e[i] = waveheight*(math.sin(wavenumber*(wavevelocity*frame*t-x*(i-0.5*n+1)))+reflectedFactor*math.sin(wavenumber*(wavevelocity*frame*t+x*(i-0.5*n+1))))
            else:
                e[i] = waveheight*(transmittedFactor*math.sin(wavenumber*(wavevelocity*frame*t-x*(i-0.5*n)*math.sqrt(d/d2))))
        ln2.set_data(xaxis2, e)
        return ln, ln2,
        
    ani = FuncAnimation(fig, update, frames=timesteps+1, init_func=init, interval=10, blit=True)

    #saves animation
    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title='MySimulation', artist='Matplotlib', comment='Animation')
    writer = FFMpegWriter(fps=15, metadata=metadata, bitrate=1000)
    """
    with writer.saving(fig, "depth.mp4", 100):
        for i in range(timesteps):
            ln.set_data(xaxis, finalh[i])
            ln2.set_data(xaxis2, e)
            writer.grab_frame()
    """
    ani.save("depthCorrect.mp4",writer=writer)
    print("reflected factor: " + str(reflectedFactor))
    print("transmitted factor: " + str(transmittedFactor))
    plt.show()

    
