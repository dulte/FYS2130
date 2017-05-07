# -*- coding: utf-8 -*-
"""
Created on Sun May 07 21:57:55 2017

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


    def calcTotalEnergy(self):
        KE = np.zeros((self.N,self.T))
        """Special cases for the end times"""
        KE[:,0] = 0.5*self.ms*(((self.y[:,0] - self.y[:,1])/(self.dt))**2)
        KE[:,-1] = 0.5*self.ms*(((self.y[:,-1] - self.y[:,-2])/(self.dt))**2)
        """KE with symmetric difference quotient"""
        for t in xrange(1,T-1):
            v = (self.y[:,t+1] - self.y[:,t-1])/(2*self.dt)
            KE[:,t] = 0.5*self.ms*(v**2)


        U = np.zeros((self.N,self.T))
        U = 0.5*self.k*(self.y[1:,:]-self.y[:-1,:])**2

        plt.plot(np.arange(self.T)*dt,np.sum(KE,axis = 0))
        plt.plot(np.arange(self.T)*dt,np.sum(U,axis = 0))
        plt.plot(np.arange(self.T)*dt,np.sum(KE,axis = 0) + np.sum(U,axis = 0))
        plt.title("Energy", fontsize = 25)
        plt.xlabel("Time [s]", fontsize = 25)
        plt.ylabel("Energy [J]", fontsize = 25)
        plt.legend(["Kinetic","Potential","Total"])
        plt.show()
        

        
if __name__ == "__main__":
    m = 0.02 #kg
    k = 10. #kg/s^2

    dt = 1*np.sqrt(m/k)

    dx = 1

    N = 200

    T = 300
    edges = "reflective"


    i = np.arange(N)

    y0 = np.sin(7*np.pi*i/float(N-1))
    y_m = np.copy(y0)

    wave = waveSolver(edges,m,k,T,N,dt,dx)

    wave.solve(y0,y_m)

    wave.calcTotalEnergy()