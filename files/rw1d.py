# 1D Random Walk
import numpy as np
import random as rd
import matplotlib.pyplot as plt

N = 100
x = np.zeros(N+1)
for n in range(N):
    c = rd.randrange(0,2)
    x[n+1] = x[n] + 1 if (c == 0) else x[n] - 1

fig = plt.figure()
plt.plot(x, '-o')
plt.show()
