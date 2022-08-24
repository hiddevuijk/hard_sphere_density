import numpy as np
import matplotlib.pyplot as plt
from sys import exit

i = 5
rho = np.loadtxt("rhoz{:d}.dat".format(i))
x = rho[:,0]
y = rho[:,1]
plt.scatter(rho[:,0], rho[:,1])


dx = x[1] - x[0]
print( dx * sum(y))

plt.xlabel(r"$z$", fontsize=15)
plt.ylabel(r"$\rho(z)$", fontsize=15,rotation=0)
plt.xlim([-3,3])
plt.ylim([0,0.8])
plt.show()
