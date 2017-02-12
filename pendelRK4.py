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



# def a(t,x,v):
#     k = 10
#     b = 0.04
#
#     L = 0
#     m = .1
#     omega = np.sqrt(k/m)
#     omega_trykt = 0.9*omega
#     V_0 = 0
#     F_trykk = V_0*np.sin(omega_trykt*t)
#     # b = 0.1
#     F = k*(L-x) - b*v
#     return (F+F_trykk)/m
#
# def a_trykt(t,x,v):
#     k = 10
#     b = 0.04
#
#     L = 0
#     m = .1
#     omega = np.sqrt(k/m)
#     omega_trykt = 1*omega
#     V_0 = 0.1
#     F_trykk = V_0*np.cos(omega_trykt*t)
#     # b = 0.1
#     F = k*(L-x) - b*v
#     if (t<25):
#         return (F+F_trykk)/m
#     else:
#         return F/m
#
#
#
# def analytisk(t):
#     k = 10
#     m = .1
#     b = .1
#     omega = np.sqrt(k/m)
#     gamma = b/(2*m)
#     omega2 = np.sqrt(omega**2-gamma**2)
#     A = .1
#     phi = -.1
#
#
#     return np.exp(-gamma*t)*A*np.cos(omega2*t+phi)
#
#
#
# T = 40
# dt = 1e-3
# x0 = np.array([0,.1])
# v0 = np.array([0,0])
#
# b = 1.5
#
# t1,x1,v1 = RK4(T,dt,x0,v0,a_trykt)
#
# oneOverE = np.max(x1[:,1])/np.exp(1)
#
# plt.plot(t1,x1[:,1])
# plt.plot(t1,np.ones(len(t1))*oneOverE)
# # plt.plot(t1,analytisk(t1))
# # plt.legend(["Numerical", "Analytic"])
# plt.xlabel("Time[s]")
# plt.ylabel("Height[m]")
# plt.title("Spring Pendulum with Force")
# plt.show()

# b = 1.5
#
# t1,x1,v1 = RK4(T,dt,x0,v0,a)
#
# b = 2
# t2,x2,v2 = RK4(T,dt,x0,v0,a)
#
# b = 4
#
# t3,x3,v3 = RK4(T,dt,x0,v0,a)
#
#
# plt.title("Behavior for different damping")
# plt.xlabel("Time[s]")
# plt.ylabel("Height[m]")
# plt.plot(t1,x1[:,1])
# plt.plot(t2,x2[:,1])
# plt.plot(t3,x3[:,1])
# plt.legend(["b = 1.5", "b = 2", "b = 4"])
# plt.show()


def a(t,x,v):
    k = 10
    b = 0.04
    L = 0
    m = .1
    # omega = 12
    V_0 = 0.1
    F_trykk = V_0*np.cos(omega*t)

    F = k*(L-x) - b*v
    # print F/m + np.array([0,-9.81])
    return (F+F_trykk)/m #+ np.array([0,-9.81])



omegas = np.linspace(5,15,200)
Emax = np.zeros_like(omegas)

T = 4
dt = 1e-2
x0 = np.array([0,.1])
v0 = np.array([0,0])

for i in range(0,len(omegas)):
    omega = omegas[i]

    t,x,v = RK4(T,dt,x0,v0,a)

    Emax[i] = np.max(0.5*0.04*np.abs(x)**2 + 0.5*0.1*np.abs(v)**2)

plt.plot(omegas/(2*np.pi),Emax)
plt.show()
