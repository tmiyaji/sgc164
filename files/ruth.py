import numpy as np
import matplotlib.pyplot as plt

def Ruth(t, q, p, dT, dU, h):
    Q1 = q + 7 * h * dT(p) / 24
    P1 = p - 2 * h * dU(Q1) / 3
    Q1 = Q1 + 0.75 * h * dT(P1)
    P1 = P1 + 2 * h * dU(Q1) / 3
    Q1 = Q1 - h * dT(P1) / 24
    P1 = P1 - h * dU(Q1)
    return [Q1, P1]

def dT(p):
    return p

def dU(q):
    r = np.hypot(q[0],q[1])
    return q / r**3

q0 = np.array([1.0, 1.0])
p0 = np.array([0.0, 1.0])
a = 0.0 # 初期時刻
h = 0.01
N = 10000
b = a + N * h

Q,P = [q0], [p0]
q,p = q0, p0
for n in range(1,N+1):
    t = a + n * h
    q,p = Ruth(t, q, p, dT, dU, h)
    Q.append(q)
    P.append(p)

QQ = np.array(Q)
PP = np.array(P)
fig = plt.figure(figsize=(5,5))
plt.plot(QQ[:,0], QQ[:,1], '-')
plt.show()
