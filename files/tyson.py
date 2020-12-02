import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def tyson(t, x, a, b, e):
    return [(x[0]*(1-x[0]) - b*x[1]*(x[0]-a)/(x[0]+a))/e, x[0] - x[1]]

a, b, e=0.01, 1.0, 0.1
t0, t1 = 0.0, 100.0
x0 = [0.1, 0.4]
s0 = solve_ivp(tyson, [t0, t1], x0, args=([a,b,e]), method='DOP853', dense_output=True)

def f(x, y):
    return (x*(1-x) - b*y*(x-a)/(x+a))/e
def g(x,y):
    return x - y

N = 20000
T = np.linspace(t0, t1, N)
sol = s0.sol(T)

fig = plt.figure()
plt.xlim(0.,0.8)
plt.ylim(0.,0.5)
Nx, Ny = 17, 17
X, Y = np.meshgrid(np.linspace(0.0, 0.8, Nx), np.linspace(0.0, 0.5, Ny))
U, V = f(X, Y), g(X, Y)
U, V = U/np.hypot(U, V), V/np.hypot(U, V)
plt.quiver(X, Y, U, V, angles='xy', color='gray')
plt.plot([0, 0.8], [0,0.8], '-')
X2 = np.linspace(0,0.8,1024)
plt.plot(X2, (X2+a)*X2*(1-X2)/(b*(X2-a)), '-')
plt.plot(sol.T[-N//10:,0], sol.T[-N//10:,1], '-')
plt.show()
