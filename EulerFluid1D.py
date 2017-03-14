import numpy as np
import matplotlib.pyplot as plt

dx = 0.01
dt = 0.01
I = int(1/dx)
N = int(1/dt)
x = np.linspace(0, 1, I-1)
t = [1.0, 2.0, 3.0]
gamma = 1.4

def rho_0(x):
    if x<= 0.5:
        return 1.0
    if x> 0.5:
        return 0.125

def P_0(x):
    if x<= 0.5:
        return 1.0
    if x> 0.5:
        return 0.1

def e_0(x):
    P_0(x)/(gamma-1)

def MacCormack():
    U1 = np.zeros([I, N])
    for i in range(I):
        U1[i][0] = rho_0(i*dx)

    U2 = np.zeros([I, N])

    U3 = np.zeros([I, N])
    for i in range(I):
        U3[i][0] = e_0(i*dx)


    F1 = np.zeros([I, N])

    F2 = np.zeros([I, N])
    for i in range(I):
        F2[i][0] = P_0(i*dx)

    F3 = np.zeros([I, N])

    
    for n in range(N-1):
        
        U1s = np.zeros(I)
        U2s = np.zeros(I)
        U3s = np.zeros(I)
        F1s = np.zeros(I)
        F2s = np.zeros(I)
        F3s = np.zeros(I)
        
        for i in range(I-1):

            U1s[i] = U1[i][n] - dt/dx*(F1[i+1][n] - F1[i][n]) 
            U2s[i] = U2[i][n] - dt/dx*(F2[i+1][n] - F2[i][n]) 
            U3s[i] = U3[i][n] - dt/dx*(F3[i+1][n] - F3[i][n]) 

        for i in range(I-1):

            F1s[i] = U2s[i]
            F2s[i] = (U2s[i])**2/U1s[i] + (gamma-1)*(U3s[i] - 0.5*(U2s[i])**2/U1s[i])
            F3s[i] = U2s[i]/U1s[i]*(U3s[i] + (gamma-1)*(U3s[i] - 0.5*(U2s[i])**2/U1s[i]))
            
        for i in range(I-1):

            U1[i][n+1] = 0.5*(U1[i][n] + U1s[i] - dt/dx*(F1s[i] - F1s[i-1]))
            U2[i][n+1] = 0.5*(U2[i][n] + U2s[i] - dt/dx*(F2s[i] - F2s[i-1]))
            U3[i][n+1] = 0.5*(U3[i][n] + U3s[i] - dt/dx*(F3s[i] - F3s[i-1]))
            
        for i in range(I-1):

            F1[i][n+1] = U2[i][n+1]
            F2[i][n+1] = (U2[i][n+1])**2/U1[i][n+1] + (gamma-1)*(U3[i][n+1] - 0.5*(U2[i][n+1])**2/U1[i][n+1])
            F3[i][n+1] = U2[i][n+1]/U1[i][n+1]*(U3[i][n+1] + (gamma-1)*(U3[i][n+1] - 0.5*(U2[i][n+1])**2/U1[i][n+1]))


    rhof = []
    for i in range(I-1):
        rhof.append(U1[i][N-1])

    Pf = []
    for i in range(I-1):
        Pf.append(F2[i][N-1]-(U2[i][N-1])**2/U1[i][N-1])

    uf = []
    for i in range(I-1):
        uf.append(U2[i][N-1]/U1[i][N-1])
    

    return rhof, Pf, uf


   
rho, P, u = MacCormack()

'''
plt.plot(x, rho, label='rho')
plt.legend()
plt.show()

plt.plot(x, P, label='P')
plt.legend()
plt.show()

plt.plot(x, u, label='u')
plt.legend()
plt.show()
'''

