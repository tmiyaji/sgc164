import numpy as np
from scipy.linalg import solve, norm
import matplotlib.pyplot as plt

# 予測子
def predict(y, h, fx, fa, w):
    n = len(w)
    A = np.zeros((n, n))
    A[:-1, :-1] = fx(y[:-1], np.array([y[-1]]))
    A[:-1, -1] = fa(y[:-1], np.array([y[-1]]))
    A[-1,:] = w
    b = np.zeros(n)
    b[-1] = 1.0
    v = solve(A, b)
    v /= norm(v)
    return y+h*v, v

# 修正子
def correct(x, h, y, v, f, fx, fa, eps=1.0e-8, maxiter=10):
    def F(x, h, y, v, f, fx, fa):
        ans = np.zeros(len(x))
        ans[:-1] = f(x[:-1],np.array([x[-1]]))
        ans[-1] = np.dot(x - y, v) - h
        return ans

    def DF(x, h, y, v, f, fx, fa):
        A = np.zeros((len(x), len(x)))
        A[:-1, :-1] = fx(x[:-1],np.array([x[-1]]))
        A[:-1, -1] = fa(x[:-1],np.array([x[-1]]))
        A[-1,:] = v
        return A
    
    z = np.copy(x)
    b = F(z, h, y, v, f, fx, fa)
    for m in range(maxiter):
        A = DF(z, h, y, v, f, fx, fa)
        z -= solve(A, b)
        b = F(z, h, y, v, f, fx, fa)
        if norm(b) < eps:
            return z
    print("not converged")

# 解の追跡バージョン0
def pathfollow00(x, a, func, dfdx, dfda, nmax=100, h=0.1, epsr=1.0e-10, epsb=1.0e-10):
    bd = []
    ndim = len(x) + len(a)
    y = np.zeros(ndim)
    y[:-1] = x
    y[-1] = a[0]
    w = np.zeros(ndim)
    w[-1] = 1.0
    x0 = y[:-1]
    a0 = np.array([y[-1]])
    # 予測子
    yt,v=predict(y, h, dfdx, dfda, w)
    bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v})

    for m in range(1, nmax+1):
        # 修正子
        y = correct(yt, h, y, v, func, dfdx, dfda, eps=epsr)
        w = v
        x0 = y[:-1]
        a0 = np.array([y[-1]])
        # 予測子
        yt, v = predict(y, h, dfdx, dfda, w)

        bd.append({'TY':"R", 'x':x0, 'a':a0, 'v':v})
    return bd


def func(x, a):
    return np.array([a[0] + x[0]**2/2 - x[0]**4])

def dfdx(x, a):
    return np.array([[x[0]-4*x[0]**3]])

def dfda(x,a):
    return np.array([1.0])

x=np.array([1.0])
a=np.array([0.5])
bd=pathfollow00(x, a, func, dfdx, dfda,nmax=200, h=-0.05)

bd_r = np.array([pt['a'][0] for pt in bd])
bd_x = np.array([pt['x'][0] for pt in bd])

fig = plt.figure()
plt.xlim(-0.2,0.5)
plt.ylim(-1.1, 1.1)
plt.xlabel(r"$\alpha$")
plt.ylabel("$x$")
plt.plot(bd_r, bd_x, 'o-')
plt.show()
