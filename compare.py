#!/usr/bin/env python

#THIS VERSION SIMULATES TSUNAMI SQUARES IN 1D

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
        a[i] += (arrayhu[i+1]-arrayhu[i])*0.1
        a[i+1] +=  (arrayhu[i+1]-arrayhu[i])*0.1*(-1)

    for i in range(n):
        arrayhu[i] += a[i]

    return arrayhu

def energy(arrh, arrv):
    global d
    global g

    offset = 0
    total = 0
    ke = 0
    pe = 0
    for i in range(len(arrh)):
        offset+=(-1)*(g*(d)**2)/2
        pe+=(arrh[i]+d)*g*((-1)*d/2+arrh[i]/2)
        ke+=(arrh[i]+d)*0.5*(arrv[i])**2
    totale = pe+ke
    print(pe-offset,ke,totale-offset)
    return totale

def move(Harr, Varr, b):
    global x
    global t
    global g
    global n

    accel = np.zeros(n)
    distmoved = np.zeros(n)
    addedheights = np.zeros(n)
    addedmomentum = np.zeros(n)
    tempmomentum = np.zeros(n)
    finalv = np.zeros(n)
    finalh = np.zeros(n)

    #calculate acceleration
    for i in range(n):
        if(i==0):
            accel[i] = (-1)*g*(Harr[i+1]-0)/(2*x)
        elif(i==n-1):
            accel[i] = (-1)*g*(0-Harr[i-1])/(2*x)
        else:
            accel[i] = (-1)*g*(Harr[i+1]-Harr[i-1])/(2*x)
    
    #calculate distance moved and new momentum
    for i in range(n):
        distmoved[i] = Varr[i]*t + 0.5*accel[i]*t**2
        tempmomentum[i] = (b[i]+Harr[i])*(Varr[i]+t*accel[i])
           
    #calculate added height and momentum
    for i in range(n):
        if(i>0 and i<n-1):
            if(distmoved[i]>0):
                addedheights[i+1] += (b[i]+Harr[i])*distmoved[i]/x
                addedheights[i] += (b[i]+Harr[i])*(1-distmoved[i]/x)
                addedmomentum[i+1] += tempmomentum[i]*distmoved[i]/x
                addedmomentum[i] += tempmomentum[i]*(1-distmoved[i]/x)
            if(distmoved[i]<=0):
                addedheights[i-1] += (b[i]+Harr[i])*(-1)*distmoved[i]/x
                addedheights[i] += (b[i]+Harr[i])*(1-(-1)*distmoved[i]/x)
                addedmomentum[i-1] += tempmomentum[i]*(-1)*distmoved[i]/x
                addedmomentum[i] += tempmomentum[i]*(1-(-1)*distmoved[i]/x)
    
    #update the heights and velocities
    for i in range(n):
        finalv[i] = addedmomentum[i]/addedheights[i]
        finalh[i] = addedheights[i]-b[i]
    finalh[0] = 0
    finalh[n-1] = 0
    finalv[0] = 0
    finalv[n-1] = 0

    return finalh, finalv
        
def gaus(mean, sd, height, where):
    
    return height*math.exp(-((where-mean)**2)/(2*sd**2))

if __name__ == "__main__":
    n = 480     #number of bins
    x = 1850.0     #distance between bins
    t = 5.0    #time step between iterations
    g = 9.81     #strength of gravity
    d = 4000.0     #depth of water
    timesteps = 50     #number of timesteps
    passes = 0     #number of smoothing passes
    tmax = x/(math.sqrt(g*d))    #maximum time step allowed

    print("maximum t: {}".format(tmax))

    if (t>tmax):
        t=tmax
        print("new t: {}".format(t))
    
    h = np.zeros(n)
    u = np.zeros(n)

    #gaussian initial condition and bathymetry
    height = 5.0
    mean = float(x*(n/2))
    sd = 5000.0
    for i in range(n):
        h[i] = gaus(mean, sd, height, i*x)

    bathymetry = np.zeros(n)
    for i in range(n):
        bathymetry[i] = d

    #calculating shallow water solver version
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

    #calculating tsunami squares version
    TSfinalh = np.zeros([timesteps+1,n])
    TSfinalh[0] = h

    TSfinalv = np.zeros([timesteps+1,n])
    TSfinalv[0] = u

    a = np.zeros(n)
    b = np.zeros(n)
    
    for i in range(timesteps):
        for j in range(0):
            TSfinalh[i] = smooth(TSfinalh[i])
            TSfinalv[i] = smooth(TSfinalv[i])
        TSfinalh[i+1], TSfinalv[i+1] = move(TSfinalh[i], TSfinalv[i], bathymetry)

    #setting up tsunami squares data
    simdata = Dataset("/home/davidgrzan/Tsunami/testing/output/test7.nc", 'r', format='NETCDF4')
    level = np.array(simdata.variables['level'])
    simdata2 = Dataset("/home/davidgrzan/Tsunami/testing/output/test8.nc", 'r', format='NETCDF4')
    level2 = np.array(simdata2.variables['level'])
    simdata3 = Dataset("/home/davidgrzan/Tsunami/testing/output/test9.nc", 'r', format='NETCDF4')
    level3 = np.array(simdata3.variables['level'])

    xaxis = np.zeros(n)
    xaxis2 = np.zeros(n)
    w = n-1
    for i in range(n):
        xaxis[i]=i*x
        xaxis2[i]=w*x
        w-=1
               
    fig, ax = plt.subplots()
    ln, = plt.plot([], [], 'k-', animated=True)
    ln2, = plt.plot([], [], 'r-', animated=True)
    ln3, = plt.plot([], [], 'g-', animated=True)
    ln4, = plt.plot([], [], 'b-', animated=True)
    ln5, = plt.plot([], [], 'b-', animated=True)
    
    def init():
        ax.set_xlim(0,(n-1)*x)
        ax.set_xlim(230*x,250*x)
        ax.set_ylim(-1,5)
        return ln, ln2, ln3, ln4, ln5,
        
    def update(frame):
        ln.set_data(xaxis, finalh[frame])
        ln5.set_data(xaxis, TSfinalh[frame])
        #ln2.set_data(xaxis2, level[frame][0])
        #ln3.set_data(xaxis2, level2[frame][0])
        #ln4.set_data(xaxis2, level3[frame][0])
        return ln, ln2, ln3, ln4, ln5,
    
    #setting up video plotter
    ani = FuncAnimation(fig, update, frames=timesteps+1, init_func=init, interval=100, blit=True)
    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title='MySimulation', artist='Matplotlib', comment='Animation')
    writer = FFMpegWriter(fps=1, metadata=metadata, bitrate=1000)
    
    with writer.saving(fig, "comparison.mp4", 100):
        for i in range(timesteps):
            ln.set_data(xaxis, finalh[i])
            #ln2.set_data(xaxis2, level[i][0])
            #ln3.set_data(xaxis2, level2[i][0])
            #ln4.set_data(xaxis2, level3[i][0])
            ln5.set_data(xaxis, TSfinalh[i])
            writer.grab_frame()
            #energy(finalh[i],finalu[i])
    
    plt.show()
