import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def f(t, x, a, b):
    dxdt = a - (b+1) * x[0] + x[1] * x[0]**2
    dydt = b * x[0] - x[1] * x[0]**2
    return np.array([dxdt, dydt])

a, b = 1.0, 2.5
t0, t1 = 0.0, 50.0
x1 = [0.1,0.1]
sol = solve_ivp(f, [t0, t1], x1, args=([a, b]), dense_output=True)

fig = plt.figure(figsize=(8,5))
plt.xlim(t0,t1)
plt.ylim(0,5)
T = np.linspace(t0,t1,512)
Z = sol.sol(T)
plt.plot(T, Z.T[:,0],'-', label="$x(t)$")
plt.plot(T, Z.T[:,1],'-', label="$y(t)$")
plt.legend()
plt.show()
