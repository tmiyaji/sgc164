import numpy as np
import random as rd
import matplotlib.pyplot as plt

def birth(beta):
    N = 10000
    a = np.zeros(N+1, dtype=int)
    a[0] = 1
    for n in range(N):
        a[n+1] = a[n] + (1 if rd.random() < beta*a[n] else 0)
    return a

# 一つのサンプルパス
rd.seed(2020202)
beta = 0.0003
a = birth(beta)
fig = plt.figure()
plt.suptitle("a sample for beta=0.0003")
plt.xlim(0, 10000)
plt.ylim(0, 40)
plt.plot(a, '.')

# 10個のサンプルパス
beta = 0.0007
a10 = np.array([birth(beta) for k in range(10)])
fig2 = plt.figure()
plt.suptitle("10 samples for beta=0.0007")
plt.xlim(0,10000)
for k in range(10):
    plt.plot(a10[k,:], 'o',ms=0.25)

# 10個のサンプルパスの平均
fig3 = plt.figure()
plt.suptitle("mean of 10 samples for beta=0.0007")
plt.xlim(0,10000)
plt.plot(np.mean(a10, axis=0), 'o', markersize=0.25)

plt.show()
