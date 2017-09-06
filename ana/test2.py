import matplotlib.pyplot as plt
import numpy as np

x,y,z, R = np.loadtxt('data.txt',usecols=(0,1,2,3), unpack = True)


fig4 = plt.figure(4)
plt.hist(x, 50)
plt.show()

fig5 = plt.figure(5)
plt.hist(y, 50)
plt.show()

fig6 = plt.figure(6)
plt.hist(z, 50)
plt.show()

fig7 = plt.figure(7)
plt.hist(R, 50)
plt.show()