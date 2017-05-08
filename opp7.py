# -*- coding: utf-8 -*-
"""
Created on Mon May 08 17:18:47 2017

@author: dulte
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import seaborn
from scipy.signal import argrelextrema


class waveSolver:

    def __init__(self, edges, m, k, T,N, dt,dx):
        
        """Checks if legal option for ending is given"""
        if edges != "open" and edges != "reflective":
            print "That is not a valid type of boundary condition!"
            exit()

        """Checks if dt fulfills the convergence conditions"""
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
        self.ms = np.ones(N)*self.m #Array of masses. Constant mass set as diffult




    def solve(self,y,y_prev,plot = True):


        M =1e10 #'infinite' mass of the endpoints

        self.y[:,0] = np.copy(y)
        
        """ Im setting y_i^- for the first timestep in the last timestep, 
        this way we dont need any special treatment of this step
        (since '-1' is the last position of an array in python)"""
        self.y[:,-1] = np.copy(y_prev)

        if self.edges == "reflective":
            self.ms[0] = float(M) #Givs the endpoint there 'infinite' mass
            self.ms[-1] = float(M)

        """This factor is predefined to save computational time"""
        self.factor = (self.dt**2/(self.ms))*self.k
                      
        
        for t in xrange(self.T-1):

            self.y[0,t+1] = (self.dt**2/float(self.ms[0]))*(self.k*self.y[1,t] - self.k*self.y[0,t] )+2*self.y[0,t] - self.y[0,t-1]
            self.y[self.N-1,t+1] = (self.dt**2/float(self.ms[self.N-1]))*(-self.k*self.y[self.N-1,t] + self.k*self.y[self.N-2,t])+2*self.y[self.N-1,t] - self.y[self.N-1,t-1]

            self.y[1:-1,t+1] =   self.factor[1:-1]*(self.y[2:,t] - 2*self.y[1:-1,t] + self.y[0:-2,t])+2*self.y[1:-1,t] - self.y[1:-1,t-1]
    
    
    """Method takes either a list of frames one wants to plot, or just a single frame.
    Set animate to true if you want to animate(is somewhat unstable, so may get errors, 
    and dont seem to work on iOS)"""
    def plotWave(self,frames,animate = False):
        if animate:
            fig = plt.figure()
            animationImages = []
            for t in range(self.T):
                animationImages.append(plt.plot(self.y[:,t], 'r'))

            ani = anim.ArtistAnimation(fig, animationImages,interval = self.dt*1000, blit = True)
            #self.fig.show()
            plt.show()
                
        try:
            for f in frames:
                if f > self.T:
                    continue
                plt.plot(self.y[:,f])
                plt.title("String at timestep %g" %f, fontsize = 25)
                plt.xlabel("Position [$\Delta x]$", fontsize = 25)
                plt.ylabel("Amplitude", fontsize = 25)
                plt.show()
        except:
            if frames > self.T:
                return
            plt.plot(self.y[:,frames])
            plt.title("String at timestep %g" %frames, fontsize = 25)
            plt.xlabel("Position [$\Delta x]$", fontsize = 25)
            plt.ylabel("Amplitude ", fontsize = 25)
            plt.show()
            
            
def intialConditions(y,N):

     for i in range(N):
         if i >= 70 and i <= 99:
             y[i] = (i-69)/30.
         elif i >= 100 and i <= 128:
             y[i] = (129- i)/30.
         else:
             y[i] = 0
              
              
if __name__ == "__main__":
    m = 0.02 #kg
    k = 10. #kg/s^2
    dt = .1*np.sqrt(m/k)
    dx = 1
    N = 200
    T = 4000
    edges = "reflective"


    y0 = np.zeros(N)
    y_m = np.zeros(N)

    intialConditions(y0,N)
    y_m =np.copy(y0)

    wave = waveSolver(edges,m,k,T,N,dt,dx)

    wave.solve(y0,y_m)

    wave.plotWave(frames = [0,100,200,500,2000,1500],animate = False)