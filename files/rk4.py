import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def RK4(t, x, f, h):
    k1 = f(t, x)
    k2 = f(t+0.5*h, x+0.5*h*k1)
    k3 = f(t+0.5*h, x+0.5*h*k2)
    k4 = f(t+h, x+h*k3)
    return x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

def func(t,x):
    return x

a, b = 0.0, 1.0
N = 10
h = (b-a)/N
t = a
x = 1.0
X = np.zeros(N+1)
X[0] = x
for n in range(N):
    x = RK4(t, x, func, h)
    X[n+1] = x
    t = a + (n + 1) * h

fig = plt.figure()
xx = np.linspace(a,b,256)
plt.plot(xx, np.exp(xx), '--')
plt.plot(np.linspace(a,b,N+1), X, 'o-')
plt.show()
