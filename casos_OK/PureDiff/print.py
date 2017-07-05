
import sys
import numpy as np
import matplotlib.pyplot as plt
from KratosMultiphysics import *    #we import the KRATOS
from KratosMultiphysics.HeatEquationApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem
import KratosMultiphysics.ExternalSolversApplication


h = np.arange(0.0,1.0+1.0/3.0,1.0/3.0)
print(h)
x = (h.shape[0])
time_step = np.arange(0.0,10.0)
# print(time_step)
t = (time_step.shape[0])
# print(t)
C = np.zeros((t,x))  # [10 x 20]
print(h.shape[0])
print((time_step.shape[0]))
# print((C.shape[0]))
# print((C.shape[1]))
z = []
# Analytical
for i in range(1,t+1):  # Analyti solution of the points [X,0.45] for different time steps
    for j in range(1,x+1):
        z = h[j-1]**2 + 0.66667**2 + 0.66667*h[j-1]*time_step[i-1]
        C[i-1][j-1] = z


 # print(C)

plt.title('Pure Diffusion')
plt.ylabel('Temperature [°C]')
plt.xlabel('Time [s]')
F = [0.4444444,0.755556,1.28889,2.04444]

# plt.grid(True)
plt.plot(h,C[9][:], 'ro',label="Analytical")
plt.plot(h,F,label="Numerical")
#plt.plot(T,C, 'ro',label="Analytical")
plt.legend(loc=2, borderaxespad=0.)
plt.show()
