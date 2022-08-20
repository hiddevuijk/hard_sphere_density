import numpy as np
import matplotlib.pyplot as plt
from sys import exit

rho = np.loadtxt("rhoz.dat")
plt.scatter(rho[:,0], rho[:,1])

#rho = np.loadtxt("rhoz1.1.dat")
#plt.scatter(rho[:,0], rho[:,1], label="1.1")
#rho = np.loadtxt("rhoz1.5.dat")
#plt.scatter(rho[:,0], rho[:,1], label="1.5")
#plt.legend(title="A")
dx = rho[1,0] - rho[0,0]
norm =  dx * sum(rho[:,1])

A = 1.5
N = 1000
x = np.linspace( rho[0,0]*1.1, rho[-1,0]*1.1, N)
y = norm * np.exp( - A * x * x) * np.sqrt( A / np.pi)
#plt.plot(x,y)

dx = x[1] - x[0]
print( dx * sum(y))

plt.xlabel(r"$z$", fontsize=15)
plt.ylabel(r"$\rho(z)$", fontsize=15,rotation=0)
plt.xlim([-3,3])
plt.ylim([0,0.8])
plt.show()
