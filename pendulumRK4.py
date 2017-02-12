import numpy as np
import matplotlib.pyplot as plt



def RK4(T,dt,x0,v0,a):

    dim = len(x0)
    timeSteps = int(float(T)/dt)
    t = np.zeros(timeSteps)
    x = np.zeros((timeSteps,dim))
    v = np.zeros((timeSteps,dim))
    x[0] = x0
    v[0] = v0


    for i in range(0,int(timeSteps)-1):

        v1 = v[i]
        a1 = a(t[i],x[i],v[i])

        v2 = v[i] + a1*dt/2.0
        a2 = a(t[i] + 0.5*dt, x[i] + v2*dt*0.5,v2)

        v3 = v[i] + a2*dt/2.0
        a3 = a(t[i] + 0.5*dt, x[i] + v3*dt*0.5,v3)

        v4 = v[i] + a3*dt
        a4 = a(t[i] + dt, x[i] + v3*dt,v4)

        aMid = 1.0/6.0*(a1 + 2*a2 + 2*a3 + a4)
        vMid = 1.0/6.0*(v1 + 2*v2 + 2*v3 + v4)


        t[i+1] = t[i] + dt
        x[i+1] = x[i] + vMid*dt
        v[i+1] = v[i] + aMid*dt

    return t,x,v



def a(t,x,v):
    k = 10

    L = 0
    m = .1
    omega = 12
    V_0 = 0.1
    F_trykk = 0#V_0*np.sin(omega*t)
    b = 0.1
    F = k*(L-x) - b*v
    return (F+F_trykk)/m


def analytisk(t):
    k = 10
    m = .1
    b = .1
    omega = np.sqrt(k/m)
    gamma = b/(2*m)
    omega2 = np.sqrt(omega**2-gamma**2)
    A = .1
    phi = -.1


    return np.exp(-gamma*t)*A*np.cos(omega2*t+phi)



T = 10
dt = 1e-3
x0 = np.array([0,.1])
v0 = np.array([0,0])

b = 1.5

t1,x1,v1 = RK4(T,dt,x0,v0,a)


plt.plot(t1,x1[:,1])
plt.plot(t1,analytisk(t1))
plt.legend(["Numerical", "Analytic"])
plt.xlabel("Time[s]")
plt.ylabel("Height[m]")
plt.title("Spring Pendulum")
plt.show()
