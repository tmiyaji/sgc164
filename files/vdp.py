import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def vdp(t, x, eps):
    return [x[1], eps * (1 - x[0]**2)*x[1] - x[0]]

t0 = 0.0
t1 = 150.0
N = 15000
tt = np.linspace(t0, t1, N)
eps = 0.1
x1 = [0.1,0.1]

s1 = solve_ivp(vdp, [t0, t1], x1, args=([eps]), t_eval=tt)
fig = plt.figure(figsize=(5,5))
plt.xlim(-4.5,4.5)
plt.ylim(-4.5,4.5)
plt.plot(s1.y[0], s1.y[1], '-')
plt.show()
