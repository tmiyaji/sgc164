import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def Euler(t, x, f, h):
    return x + h * f(t, x)

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
    x = Euler(t, x, func, h)
    X[n+1] = x
    t = a + (n + 1) * h

fig = plt.figure()
xx = np.linspace(a,b,256)
plt.plot(xx, np.exp(xx), '--')
plt.plot(np.linspace(a,b,N+1), X, 'o-')
plt.show()
