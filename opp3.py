# -*- coding: utf-8 -*-
"""
Created on Mon May 01 14:29:04 2017

@author: dulte
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim


class waveSolver:
    
    def __init__(self, edges, m, k, N, T, dt,dx):
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
        self.T = T
        self.dx = dx
        self.dt = dt
        self.validInitialConditions = False
        
        
        self.y = np.zeros((T,N))
        
        
        
    def setInitialConditions(self, y0, **kwargs):
#        if (type(kwargs.get("v0",None)) is not type(None)):
#            self.v[0,:] = kwargs.get("v0",None)
#        if (type(kwargs.get("v_m",None)) is not type(None)):
#            self.v_m = kwargs.get("v_m",None)
#        else:
#            self.v_m = np.zeros(self.N)

        if (type(kwargs.get("y_m",None)) is not type(None)):
            self.y_m = kwargs.get("y_m",None)
        else:
            self.y_m = np.zeros(self.N)
        
        self.y[0,:] = y0
        self.validInitialConditions = True
        

    
    def solve(self,animate = False):
        
        if not self.validInitialConditions:
            print "No valid initial conditions were given!"
            return
        
        M =1e10
        
        if animate:
            self.fig = plt.figure()
            self.animationImages = []
            
        
        if self.edges == "open":
            self.y[1,0] = (self.dt**2/float(self.m))*(self.k*self.y[0,1] - 2*self.k*self.y[0,0] + self.k*self.y[0,N-1]) + 2*self.y[0,0] - self.y_m[0]
            self.y[1,N-1] = (self.dt**2/float(self.m))*(self.k*self.y[0,0] - 2*self.k*self.y[0,N-1] + self.k*self.y[0,N-2]) + 2*self.y[0,N-1] - self.y_m[N-1]
            
            
        elif self.edges == "reflective":
            self.y[1,0] = (self.dt**2/float(M))*(self.k*self.y[0,1] - 2*self.k*self.y[0,0] + self.k*self.y[0,N-1]) + 2*self.y[0,0] - self.y_m[0]
            self.y[1,N-1] = (self.dt**2/float(M))*(self.k*self.y[0,0] - 2*self.k*self.y[0,N-1] + self.k*self.y[0,N-2]) + 2*self.y[0,N-1] - self.y_m[N-1]
           
            
            
        for i in range(1,N-1):
            self.y[1,i] = (self.dt**2/float(self.m))*(self.k*self.y[0,i+1] - 2*self.k*self.y[0,i] + self.k*self.y[0,i-1]) + 2*self.y[0,i] - self.y_m[i]

                
        for t in range(1,T-1):
            
            if self.edges == "open":
                self.y[t+1,0] = (self.dt**2/float(self.m))*(self.k*self.y[t,1] - 2*self.k*self.y[t,0] + self.k*self.y[t,N-1]) + 2*self.y[t,0] - self.y[t-1,0]
                self.y[t+1,N-1] = (self.dt**2/float(self.m))*(self.k*self.y[t,0] - 2*self.k*self.y[t,N-1] + self.k*self.y[t,N-2]) + 2*self.y[t,N-1] - self.y[t-1,N-1]
               
                
                
            elif self.edges == "reflective":
                self.y[t+1,0] = (self.dt**2/float(M))*(self.k*self.y[t,1] - 2*self.k*self.y[t,0] + self.k*self.y[t,N-1]) + 2*self.y[t,0] - self.y[t-1,0]
                self.y[t+1,N-1] = (self.dt**2/float(M))*(self.k*self.y[t,0] - 2*self.k*self.y[t,N-1] + self.k*self.y[t,N-2]) + 2*self.y[t,N-1] - self.y[t-1,N-1]
                
                
                
            for i in range(1,N-1):
                self.y[t+1,i] = (self.dt**2/float(self.m))*(self.k*self.y[t,i+1] - 2*self.k*self.y[t,i] + self.k*self.y[t,i-1]) + 2*self.y[t,i] - self.y[t-1,i]
    
            
            
            
            if animate:
                self.animationImages.append(plt.plot(self.y[t,:],"r"))
                
    
    def plotWave(self,frame = 0, animate = False):
        if animate and self.animationImages == None:
            print "No animation was made!"
            return
        
        if frame == 0 and not animate:
            print "No option for plotting was made!"
            return
        
        
        if animate:
            ani = anim.ArtistAnimation(self.fig, self.animationImages,interval = 50, blit = False)
            plt.show()
        
        if frame != 0:
            try:
                for f in frame:
                    plt.plot(self.y[f,:])
                    plt.title("Spring at timestep %g" %f)
                    plt.xlabel("Position")
                    plt.ylabel("Amplitude")
                    plt.show()
            
            except:
                plt.plot(self.y[frame,:])
                plt.title("Spring at timestep %g" %frame)
                plt.xlabel("Position")
                plt.ylabel("Amplitude")
                plt.show()
            
            
    

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
    
    T = 1200
    
    edges = "reflective"
    
    
    i = np.arange(N)
    
    #y0 = np.sin(7*np.pi*i/float(N-1))
    y0 = np.zeros(N)
    y_m = np.zeros(N)
    
    intialConditions(y0,N)
    prevInitialConditions(y_m,y0,N,dt)
    
    
    wave = waveSolver(edges,m,k,N,T,dt,dx)
    
    wave.setInitialConditions(y0, y_m = y_m)
    
    wave.solve(animate=True)
    wave.plotWave(frame=[0,1,2,3],animate=True)
    

   
            