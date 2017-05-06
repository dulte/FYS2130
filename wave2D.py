# -*- coding: utf-8 -*-
"""
Created on Mon May 01 14:29:04 2017

@author: dulte
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import seaborn
from scipy.signal import argrelextrema


class waveSolver:

    def __init__(self, edges, m, k, T,N, dt,dx):
        if edges != "open" and edges != "reflective":
            print "That is not a valid type of boundary condition!"
            exit()

        if dt > np.sqrt(float(m)/float(k)):
            print "m,k and dt are chosen in a way that makes dt to large!"
            exit()

        self.edges = edges
        self.m = m
        self.k = k
        self.N = N
        self.dx = dx
        self.dt = dt
        self.T = T

        self.y = np.zeros((N,T))
        self.ms = np.ones(N)*self.m







    def solve(self,y,y_prev,plot = True):


        M =1e10

        self.y[:,0] = np.copy(y)
        self.y[:,-1] = np.copy(y_prev)

        if self.edges == "reflective":
            self.ms[0] = float(M)
            self.ms[-1] = float(M)


        self.factor = (self.dt**2/(self.ms))*self.k

        for t in xrange(self.T-1):


            # self.y[0,t+1] = self.factor[0]*(self.k*self.y[1,t] - self.k*self.y[0,t] )+2*self.y[0,t] - self.y[0,t-1]
            # self.y[-1,t+1] = self.factor[-1]*(-self.k*self.y[N-1,t] + self.k*self.y[N-2,t])+2*self.y[N-1,t] - self.y[N-1,t-1]


            self.y[0,t+1] = (self.dt**2/float(self.ms[0]))*(self.k*self.y[1,t] - self.k*self.y[0,t] )+2*self.y[0,t] - self.y[0,t-1]
            self.y[N-1,t+1] = (self.dt**2/float(self.ms[N-1]))*(-self.k*self.y[N-1,t] + self.k*self.y[N-2,t])+2*self.y[N-1,t] - self.y[N-1,t-1]


            self.y[1:-1,t+1] =   self.factor[1:-1]*(self.y[2:,t] - 2*self.y[1:-1,t] + self.y[0:-2,t])+2*self.y[1:-1,t] - self.y[1:-1,t-1]






    def plotWave(self,frames,animate = False):
        if animate:
            fig = plt.figure()
            animationImages = []
            for t in range(T):
                animationImages.append(plt.plot(self.y[:,t], 'r'))

            ani = anim.ArtistAnimation(fig, animationImages,interval = self.dt*1000, blit = True)
            #self.fig.show()
            plt.show()

        try:
            for f in frames:
                if f > self.T:
                    continue
                plt.plot(self.y[:,f])
                plt.title("Spring at timestep %g" %f)
                plt.xlabel("Position")
                plt.ylabel("Amplitude")
                plt.show()
        except:
            if frames > self.T:
                return
            plt.plot(self.y[:,frames])
            plt.title("Spring at timestep %g" %frames)
            plt.xlabel("Position")
            plt.ylabel("Amplitude")
            plt.show()



    def findImpedance(self, timestep, massFactor):
        R = np.min(self.y[:,timestep])/np.max(self.y[:,0])
        analyt = (1-R)/(1+R)
        print "Numerical Impedance Difference = ", analyt#1-2*(np.min(self.y[:,timestep])/np.max(self.y[:,-1]))
        print "Analytic = ", np.sqrt(massFactor)



    def findFreq(self,index):
        y_indexed = self.y[int(index),:]
        plt.plot(np.arange(self.T)*self.dt,y_indexed)
        plt.title("Oscillation of node %g" %(index + 1))
        plt.show()
        print 1./(200.0/3.5)*np.sqrt(self.k/self.m)
        print "The frequency found analytically is %g" %(1./((self.N-1)/3.5)*np.sqrt(self.k/self.m))

        argLocalMax = argrelextrema(y_indexed,np.greater)[0]
        freqFromWaveLength = np.mean(1./(np.abs(self.dt*argLocalMax[:-1] - self.dt*argLocalMax[1:])))
        print "The frequency found from the wave length is %g" %freqFromWaveLength

        fourier = np.fft.fft(y_indexed)
        freq = np.fft.fftfreq(self.T, d = self.dt)

        print "The frequency found from the Fourier Transform is %g" %(np.abs(freq[np.argmax(np.abs(fourier))]))
        plt.plot(np.abs(freq),np.abs(fourier))
        plt.show()


    def makeDifferentMass(self,diffMassFromIndex, massFactor):
        self.ms[diffMassFromIndex:] = massFactor*self.m



    def calcTotalEnergy(self):
        KE = np.zeros((self.N,self.T))
        KE[:,0] = 0.5*self.ms*(((self.y[:,0] - self.y[:,1])/(self.dt))**2)
        KE[:,-1] = 0.5*self.ms*(((self.y[:,-1] - self.y[:,-2])/(self.dt))**2)
        for t in xrange(1,T-1):
            v = (self.y[:,t+1] - self.y[:,t-1])/(2*self.dt)
            KE[:,t] = 0.5*self.ms*(v**2)


        U = np.zeros((self.N,self.T))
        U[0,:] = 0.5*self.k*(self.y[1]-self.y[0])**2
        U[-1,:] = 0.5*self.k*(self.y[-1]-self.y[-2])**2

        U[1:-1,:] = 0.5*self.k*((self.y[2:,:]-self.y[1:-1,:])**2 +(self.y[1:-1,:]-self.y[0:-2,:])**2)/2.

        plt.plot(np.arange(self.T)*dt,np.sum(KE,axis = 0))
        plt.plot(np.arange(self.T)*dt,np.sum(U,axis = 0))
        plt.plot(np.arange(self.T)*dt,np.sum(KE,axis = 0) + np.sum(U,axis = 0))
        plt.title("Energy")
        plt.xlabel("Time [s]", fontsize = 25)
        plt.ylabel("Energy [J]", fontsize = 25)
        plt.legend(["Kinetic","Potential","Total"])
        plt.show()





def intialConditions(y,N):

    for i in range(N):

        if i >= 2 and i <= 31:
            y[i] = (i-1)/30.
        elif i >= 32 and i <= 60:
            y[i] = (61- i)/30.
        else:
            y[i] = 0
#
# def intialConditions(y,N):
#
#     for i in range(N):
#
#         if i >= 70 and i <= 99:
#             y[i] = (i-69)/30.
#         elif i >= 100 and i <= 128:
#             y[i] = (129- i)/30.
#         else:
#             y[i] = 0





if __name__ == "__main__":
    m = 0.02 #kg
    k = 10. #kg/s^2

    dt = 1*np.sqrt(m/k)
    #dt = 0.01
    dx = .001

    N = 200

    #T = 1500 #To find impedance
    T = 5100
    #T = 150
    edges = "reflective"


    i = np.arange(N)

    y0 = np.sin(7*np.pi*i/float(N-1))
    #y0 = np.zeros(N)
    #y_m = np.zeros(N)

    #intialConditions(y0,N)

    y_m =np.copy(y0)

    #y_m[1:61] = np.copy(y0[2:62])





    wave = waveSolver(edges,m,k,T,N,dt,dx)
    #wave.makeDifferentMass(100,3)



    wave.solve(y0,y_m)
    wave.plotWave(frames = [140],animate = True)

    wave.calcTotalEnergy()
    wave.findImpedance(135,3)
    wave.findFreq(99)
