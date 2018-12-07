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
        a[i] += (arrayhu[i+1]-arrayhu[i])*0.15*0.5
        a[i+1] +=  (arrayhu[i+1]-arrayhu[i])*0.15*0.5*(-1)

    for i in range(n):
        arrayhu[i] += a[i]

    return arrayhu

if __name__ == "__main__":
    n = 100     #number of bins
    x = 5585.0     #distance between bins
    t = 10    #time step between iterations
    g = 9.81     #strength of gravity
    d = 1000.0     #depth of water
    timesteps = 100     #number of iterations
    passes = 1    #number of smoothing passes
    tmax = x/(math.sqrt(g*d))     #maximum time step allowed

    print("maximum t: {}".format(tmax))

    if (t>tmax):
        t=tmax
        print("new t: {}".format(t))
    
    h = np.zeros(n)
    u = np.zeros(n)

    for i in range(n):
        if i<(n/2+5) and i>(n/2-5):
            h[i] = 0
        if (i<=(n/2-5) and i>(n/2-10)) or (i>=(n/2+5) and i<(n/2+10)):
            h[i] = 0
        if i==75:
            h[i] = 1

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


    simdata = Dataset("/home/davidgrzan/Tsunami/TestFiles/output1/testwaveFLAT.nc", 'r', format='NETCDF4')
    times = np.array(simdata.variables['time'])
    lats = np.array(simdata.variables['latitude'])
    lons = np.array(simdata.variables['longitude'])
    heightt = np.array(simdata.variables['height'])
    level = np.array(simdata.variables['level'])
    alt = simdata.variables['altitude']

    #print(times[0],times[1],times[0]-times[1],len(times),times[100])
    #print(lons.min(),lons.max(),lats.min(),lats.max(),len(lats),len(lons),lats[0]-lats[1])
    #print(level.max())
                     
    
            
    xaxis = np.zeros(n)
    xaxis2 = np.zeros(n)
    w = n-1
    for i in range(n):
        xaxis[i]=i*x
        xaxis2[i]=w*x
        w-=1
               
    fig, ax = plt.subplots()
    ln, = plt.plot([], [], 'b-', animated=True)
    ln2, = plt.plot([], [], 'r-', animated=True)
    
    def init():
        ax.set_xlim(0,(n-1)*x)
        ax.set_ylim(-0.5,0.5)
        return ln, ln2, 
        
    def update(frame):
        ln.set_data(xaxis, finalh[frame])
        ln2.set_data(xaxis2, level[frame][50])
        return ln, ln2,
        
    ani = FuncAnimation(fig, update, frames=timesteps+1, init_func=init, interval=100, blit=True)

    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title='MySimulation', artist='Matplotlib', comment='Animation')
    writer = FFMpegWriter(fps=15, metadata=metadata, bitrate=1000)

    with writer.saving(fig, "comparison.mp4", 100):
        for i in range(timesteps):
            ln.set_data(xaxis, finalh[i])
            ln2.set_data(xaxis2, level[i][50])
            writer.grab_frame()

    plt.show()
