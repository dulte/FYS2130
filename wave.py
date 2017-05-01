# -*- coding: utf-8 -*-
"""
Created on Mon May 01 14:29:04 2017

@author: dulte
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim


class waveSolver:
    
    def __init__(self, edges, m, k, N, dt,dx):
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
 
        
        
        self.y = np.zeros(N)
        self.y_prev = np.zeros(N)
        self.y_next = np.zeros(N)
        
        
    
        

    
    def solve(self,y,y_prev,T,plot = True, animate = False):
        
   
        M =1e10
        
        if animate:
            self.fig = plt.figure()
            self.animationImages = []
        
        self.y = np.copy(y)
        self.y_prev = np.copy(y_prev)
        
        self.ms = np.ones(N)*self.m
        self.ms[100:] = 3*self.m
        if self.edges == "reflective":
            self.ms[0] = M
            self.ms[-1] = M
        
        
        
        self.totalEnergy = np.zeros(T)
        
        self.factor = (self.dt**2/(self.ms))*self.k
        
        for t in range(T-1):
            
            
            self.y_next[0] = (self.dt**2/float(self.ms[0]))*(self.k*self.y[1] - 2*self.k*self.y[0] + self.k*self.y[self.N-1])+2*self.y[0] - self.y_prev[0]
            self.y_next[N-1] = (self.dt**2/float(self.ms[N-1]))*(self.k*self.y[0] - 2*self.k*self.y[N-1] + self.k*self.y[N-2])+2*self.y[N-1] - self.y_prev[N-1]
            
                
            for i in range(1,self.N-1):

                self.y_next[i] = self.factor[i]*(self.y[i+1] - 2*self.y[i] + self.y[i-1])+2*self.y[i] - self.y_prev[i]
            
            
            self.totalEnergy[t] = self.calcTotalEnergy(self.y,self.y_prev)
            
            self.y_prev = np.copy(self.y)
            self.y = np.copy(self.y_next)
            
            
            
            if animate:
                self.animationImages.append(plt.plot(self.y,"r"))
        
        
        
                
        if animate:
            ani = anim.ArtistAnimation(self.fig, self.animationImages,interval = 10, blit = False)
            plt.show()
        
        if plot:
            plt.plot(self.y)
            plt.title("Spring at timestep %g" %T)
            plt.xlabel("Position")
            plt.ylabel("Amplitude")
            plt.show()
        
        plt.plot(self.totalEnergy)
        plt.show()
            
    
    def calcTotalEnergy(self,y,y_prev):
        
        v = (y-y_prev)/dt
        KE = 0.5*m*v**2
        
        U = np.zeros(N)
        U[0] = 0.5*self.k*(y[1]-y[0])**2
        U[-1] = 0.5*self.k*(y[-1]-y[-2])**2
        for i in range(1,N-1):
            U[i] = 0.5*self.k*((y[i+1]-y[i])**2 +(y[i]-y[i-1])**2)
            
        
        return np.sum(KE) + np.sum(U)
        
        
        
    
            
            
    

def intialConditions(y,N):
    
    for i in range(N):
        
        if i >= 1 and i <= 30:
            y[i] = (i)/30.
        elif i >= 31 and i <= 59:
            y[i] = (60 - i)/30.
        else:
            y[i] = 0
             
def prevInitialConditions(y,y0,N,dt):
    for i in range(N):
        
        if i >= 1 and i <= 30:
            y[i] = y0[i] + dt*1/30.
        elif i >= 31 and i <= 59:
            y[i] = y0[i] + dt*(-1)/30.
        else:
            y[i] = 0

 

if __name__ == "__main__":
    m = 0.02 #kg
    k = 10. #kg/s^2
    
    dt = 0.1*np.sqrt(m/k)
    
    dx = 0.0001
    
    N = 200
    
    T = 2*1200
    
    edges = "reflective"
    
    
    i = np.arange(N)
    
    y0 = np.sin(7*np.pi*i/float(N-1))
    #y0 = np.zeros(N)
    y_m = np.zeros(N)
    
    #intialConditions(y0,N)
    #prevInitialConditions(y_m,y0,N,dt)
    y_m = np.copy(y0)
    
    wave = waveSolver(edges,m,k,N,dt,dx)
    
  
    
    wave.solve(y0,y_m,T,animate=True)

    

   
            
