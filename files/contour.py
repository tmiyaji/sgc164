import numpy as np
import matplotlib.pyplot as plt

def E(x, y):
    return 0.5 * y**2 - np.cos(x)

x = np.linspace(-2*np.pi,2*np.pi,512)
y = np.linspace(-2.5,2.5,512)
X, Y = np.meshgrid(x, y)
Z = E(X, Y)
fig = plt.figure()
plt.contour(X, Y, Z, levels=10)
plt.show()
